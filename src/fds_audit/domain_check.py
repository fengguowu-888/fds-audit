"""
domain_check.py — FDS 输入文件几何硬校验器
==========================================
为什么需要它
------------
FDS 对**被网格裁掉的实体**、**贴在实心格子上的 VENT** 都是**静默处理**：不报 warning、
不在 .out 里留痕，照跑照出结果。只有事后再去核几何才能发现。论文4 的 C2/C3 两套算例
就是这么把配电柜、人孔、地板、端墙整批丢掉的（详见随附论文 paper4 的 C2/C3 算例）。

本模块把"读 .fds → 重建格子 → 断言"做成一次可复用的检查，两条用法：

    python domain_check.py  <路径>            # 目录则递归找 *.fds
    from domain_check import check_file        # 生成器里调用，不通过就不落盘

检查项（每项都有明确判据，不靠 FDS 自报）
----------------------------------------
E1 containment      OBST/VENT 落在 &MESH 外的体积分数；0% = 被静默丢弃（错误）
E2 clipping         部分落在域外 = 被裁（警告，附裁掉百分比）
G1 grid_conformance 每个面到最近网格平面的偏差（以格子为单位）；>0.5 说明被 FDS 挪了位置
V1 vent_gas_path    VENT 面片两侧有没有气相格子；两侧全实心 = 该开口失效（错误）
V2 vent_area        开口实测生效面积 vs 标称面积；部分被堵（<95%）报警
O1 overlap          两个 OBST 完全重合 / 一个完全埋在另一个里面
O2 occlusion        某 OBST 覆盖的格子全部已被其它构件占用（≥99%）→ 对流动零贡献
D1 devc_in_gas      &DEVC XYZ 落在实心格子里（读数无意义）
B1 boundary_audit   网格边界上有多少气相格子、其中多少被 VENT 覆盖；
                    剩下的就是**本该是混凝土、实际退化成 INERT 的面**
R1 resolution       D*/dx（据 HRRPUA×火源面积反推 HRR）
S1 snapping_report  被 FDS 静默吸附的几何面清单（declared → snapped 坐标，离散化程度）

判据出处
--------
- FDS 把 OBST 转格子：格子**中心**落在 OBST 内部则该格子为实心。
- VENT 必须是某个格子面；落在实心格子之间的 VENT 不产生任何流动。

多网格
------
- 平铺的多个 &MESH（MPI 分解常见形态）会合并成单一逻辑网格精确校验；
  非均匀 / 有间隙 / 有重叠的多网格明确降级告警（MESH:WARN），不假装支持。
"""

from __future__ import annotations

import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ──────────────────────────────────────────────────────────────────
# 解析
# ──────────────────────────────────────────────────────────────────

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"
XB_RE = re.compile(r"\bXB\s*=\s*(" + r"\s*,\s*".join([_NUM] * 6) + r")")
IJK_RE = re.compile(r"\bIJK\s*=\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)")
XYZ_RE = re.compile(r"\bXYZ\s*=\s*(" + r"\s*,\s*".join([_NUM] * 3) + r")")
STR_RE = re.compile(r"(\w+)\s*=\s*'([^']*)'")
BOOL_RE = re.compile(r"\bTHICKEN\s*=\s*\.?(TRUE|FALSE|T|F)\.?", re.I)

# D1「落在实心格子里」只对**气相场点量**成立。墙面/表面量（壁温、热流）本来就要贴
# 在实心表面，非空间量（TIME/CONTROL/TIMER）根本不读气相，落实心不构成缺陷。
_WALL_Q = {
    "WALL TEMPERATURE", "BACK WALL TEMPERATURE", "INSIDE WALL TEMPERATURE",
    "INCIDENT HEAT FLUX", "RADIATIVE HEAT FLUX", "NET HEAT FLUX",
    "CONVECTIVE HEAT FLUX", "TOTAL HEAT FLUX", "GAUGE HEAT FLUX",
    "WALL PRESSURE", "NORMAL VELOCITY", "RADIOMETER", "GAUGE TEMPERATURE",
    "THERMOCOUPLE",
}
_NONSPATIAL_Q = {"TIME", "CONTROL VALUE", "CONTROL", "TIMER"}


def _f(s: str) -> float:
    return float(s.replace("d", "e").replace("D", "E"))


def _six(s: str):
    return tuple(_f(v) for v in s.split(","))


def _nint(x: float) -> int:
    """Fortran NINT 四舍五入：.5 向远离零取整（FDS 半格吸附的约定），带浮点容差。

    网格面坐标 = origin + n*d，n 恰为半整数（面落在格子中心）时，FDS 用 NINT 向
    远离零取整；Python 的 round(.5)=偶数 与浮点误差（1.5 实为 1.4999…）都会把
    半格吸错方向，导致「有效开口」被误判成「两侧全实心」。
    """
    a = abs(x)
    f = a - math.floor(a)
    if abs(f - 0.5) < 1e-9:
        n = math.floor(a) + 1
    else:
        n = int(math.floor(a + 0.5))
    return n if x >= 0 else -n


def parse_records(txt: str):
    """按 FDS namelist 规则切分记录：'&NAME ...' 起、'/' 止为一个记录。
    FDS 会忽略记录之外的自由文本（纯文本说明行、行首 '#' 注释），
    行内 '!' 是行注释。这比单纯按 '/' 切分更能复现 FDS 的真实读取行为。"""
    recs = []
    cur = None
    for ln in txt.splitlines():
        body = ln.split("!")[0]
        if cur is None:
            if body.lstrip().startswith("&"):
                cur = body
                if "/" in cur:
                    seg = cur.split("/", 1)[0].strip()
                    if seg:
                        recs.append(seg)
                    cur = None
        else:
            cur += "\n" + body
            if "/" in cur:
                seg = cur.split("/", 1)[0].strip()
                if seg:
                    recs.append(seg)
                cur = None
    return recs


def _parse_mult(txt: str) -> dict:
    """解析 &MULT 网格复制定义，返回 {MULT_ID: {dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi}}。

    FDS 的 &MULT 用 DX/DY/DZ 表示相邻副本的**物理偏移**（米），I/J/K_LOWER/UPPER
    表示副本索引范围（默认 0）。MESH 带 MULT_ID 时，按此规则平铺复制。
    """
    out = {}
    for rec in parse_records(txt):
        head = rec.split(None, 1)[0].upper() if rec.split() else ""
        if head != "&MULT":
            continue
        s = dict(STR_RE.findall(rec))
        mid = s.get("ID")
        if not mid:
            continue

        def num(name, default=0.0):
            mm = re.search(r"\b" + name + r"\s*=\s*(" + _NUM + r")", rec)
            return _f(mm.group(1)) if mm else default

        def inum(name, default=0):
            mm = re.search(r"\b" + name + r"\s*=\s*(-?\d+)", rec)
            return int(mm.group(1)) if mm else default

        out[mid] = dict(dx=num("DX"), dy=num("DY"), dz=num("DZ"),
                        ilo=inum("I_LOWER"), ihi=inum("I_UPPER"),
                        jlo=inum("J_LOWER"), jhi=inum("J_UPPER"),
                        klo=inum("K_LOWER"), khi=inum("K_UPPER"))
    return out


def _point_in_xb(pt, xb) -> bool:
    """点是否落在网格 XB 内（含边界）。"""
    x, y, z = pt
    return (xb[0] - 1e-9 <= x <= xb[1] + 1e-9
            and xb[2] - 1e-9 <= y <= xb[3] + 1e-9
            and xb[4] - 1e-9 <= z <= xb[5] + 1e-9)


@dataclass
class Box:
    kind: str            # 'OBST' | 'VENT'
    xb: tuple
    surf: str = ""
    thicken: bool = False
    permit_hole: bool = False  # OBST 显式 PERMIT_HOLE=F → 有意嵌套实体（免 O2 误报）
    ident: str = ""
    seq: int = 0          # 构件在文件里的出现序号（无 ID 时用这个标识）

    @property
    def name(self) -> str:
        return self.ident or f"{self.kind}#{self.seq}"

    def brief(self) -> str:
        return f"{self.name}(surf={self.surf or 'OPEN'})"

    def is_plane(self):
        """返回退化轴索引（0/1/2），非退化返回 None。"""
        deg = [i for i in range(3) if abs(self.xb[2 * i + 1] - self.xb[2 * i]) < 1e-9]
        return deg[0] if len(deg) == 1 else (tuple(deg) if deg else None)


@dataclass
class FdsModel:
    path: Path
    chid: str = ""
    ijk: tuple = ()
    mesh_xb: tuple = ()
    meshes: list = field(default_factory=list)   # [(ijk, xb), ...] 所有 &MESH（多网格）
    transformed: bool = False                    # 含 TRNX/Y/Z_ID 变换网格（Cartesian 近似无效）
    obst: list = field(default_factory=list)
    vents: list = field(default_factory=list)
    holes: list = field(default_factory=list)   # &HOLE 的 XB：从 OBST 里抠掉气相区
    burning_surfs: set = field(default_factory=set)  # 燃烧/释气面 SURF ID（HRRPUA/MLRPUA/MASS_FLUX）
    material_surfs: set = field(default_factory=set)  # 真实材质面 SURF ID（带 MATL_ID 或 THICKNESS）
    devc: list = field(default_factory=list)   # (id, quantity, xyz)
    hrrpua_sum_area: float = 0.0               # Σ HRRPUA×面积（粗估总 HRR）
    t_end: float = 0.0

    @property
    def dx(self):
        return (self.mesh_xb[1] - self.mesh_xb[0]) / self.ijk[0]

    @property
    def dy(self):
        return (self.mesh_xb[3] - self.mesh_xb[2]) / self.ijk[1]

    @property
    def dz(self):
        return (self.mesh_xb[5] - self.mesh_xb[4]) / self.ijk[2]

    @property
    def cells(self):
        return self.ijk[0] * self.ijk[1] * self.ijk[2]


def parse_fds(path: Path) -> FdsModel:
    """从磁盘文件解析。"""
    return parse_text(path.read_text(encoding="utf-8", errors="replace"), path)


def parse_text(txt: str, path: Path = None) -> FdsModel:
    """从字符串解析 —— 生成器可以在落盘之前先自检。"""
    path = path or Path("<string>")
    m = FdsModel(path=path)
    surf_hrrpua = {}

    # 先收集 SURF 属性：HRRPUA 供火源识别；燃烧面(HRRPUA/MLRPUA/MASS_FLUX)供 O2 判据
    for rec in parse_records(txt):
        if rec.upper().startswith("&SURF"):
            s = dict(STR_RE.findall(rec))
            if "ID" not in s:
                continue
            sid = s["ID"]
            h = re.search(r"HRRPUA\s*=\s*(" + _NUM + r")", rec)
            if h:
                surf_hrrpua[sid] = _f(h.group(1))
            if re.search(r"\b(?:HRRPUA|MLRPUA|MASS_FLUX)\s*=", rec):
                m.burning_surfs.add(sid)
            if re.search(r"\b(?:MATL_ID|THICKNESS)\s*=", rec):
                m.material_surfs.add(sid)

    mult = _parse_mult(txt)
    raw_meshes = []   # (ijk, xb, mult_id, transformed)
    for rec in parse_records(txt):
        head = rec.split(None, 1)[0].upper() if rec.split() else ""
        if head == "&MESH":
            mi = IJK_RE.search(rec)
            mx = XB_RE.search(rec)
            if mi and mx:
                s = dict(STR_RE.findall(rec))
                up = rec.upper()
                trn = ("TRNX_ID" in up) or ("TRNY_ID" in up) or ("TRNZ_ID" in up)
                raw_meshes.append((tuple(int(g) for g in mi.groups()),
                                   _six(mx.group(1)), s.get("MULT_ID", ""), trn))
        elif head == "&HEAD":
            s = dict(STR_RE.findall(rec))
            m.chid = s.get("CHID", "")
        elif head == "&TIME":
            t = re.search(r"T_END\s*=\s*(" + _NUM + r")", rec)
            if t:
                m.t_end = _f(t.group(1))

    # 展开 MULT 复制的网格（MULT_ID → 平铺副本集）；变换网格标记为不支持
    for ijk, xb, mid, trn in raw_meshes:
        if trn:
            m.transformed = True
            m.meshes.append((ijk, xb))
            continue
        p = mult.get(mid) if mid else None
        if p is None:
            m.meshes.append((ijk, xb))
            continue
        for i in range(p["ilo"], p["ihi"] + 1):
            for j in range(p["jlo"], p["jhi"] + 1):
                for k in range(p["klo"], p["khi"] + 1):
                    m.meshes.append((ijk, (xb[0] + i * p["dx"], xb[1] + i * p["dx"],
                                           xb[2] + j * p["dy"], xb[3] + j * p["dy"],
                                           xb[4] + k * p["dz"], xb[5] + k * p["dz"])))
    if m.meshes:
        m.ijk, m.mesh_xb = m.meshes[0]

    seq = 0
    for rec in parse_records(txt):
        head = rec.split(None, 1)[0].upper() if rec.split() else ""
        if head not in ("&OBST", "&VENT", "&HOLE", "&DEVC"):
            continue
        mx = XB_RE.search(rec)
        if head == "&DEVC":
            xy = XYZ_RE.search(rec)
            s = dict(STR_RE.findall(rec))
            if xy:
                m.devc.append((s.get("ID", ""), s.get("QUANTITY", ""),
                               tuple(_f(v) for v in xy.group(1).split(","))))
            continue
        if not mx:
            continue
        xb = _six(mx.group(1))
        if head == "&HOLE":
            m.holes.append(xb)
            continue
        seq += 1
        s = dict(STR_RE.findall(rec))
        surf = s.get("SURF_ID", "") or s.get("SURF_IDS", "")
        bx = Box(kind=head[1:], xb=xb, surf=surf, seq=seq,
                 thicken=bool(BOOL_RE.search(rec)),
                 permit_hole=bool(re.search(r"\bPERMIT_HOLE\s*=\s*\.?(?:FALSE|F)\b", rec, re.I)),
                 ident=s.get("ID", ""))
        if head == "&OBST":
            m.obst.append(bx)
        else:
            m.vents.append(bx)
            if surf in surf_hrrpua:
                deg = [i for i in range(3) if abs(xb[2 * i + 1] - xb[2 * i]) < 1e-9]
                ext = [abs(xb[2 * i + 1] - xb[2 * i]) for i in range(3) if i not in deg]
                if len(ext) == 2:
                    m.hrrpua_sum_area += surf_hrrpua[surf] * ext[0] * ext[1]

    return m


# ──────────────────────────────────────────────────────────────────
# 格子占用
# ──────────────────────────────────────────────────────────────────

class Grid:
    """按 FDS 规则重建实心格子：格子中心落在任一 OBST 内 → 实心。"""

    def __init__(self, model: FdsModel):
        self.m = model
        self.nx, self.ny, self.nz = model.ijk
        self.dx, self.dy, self.dz = model.dx, model.dy, model.dz
        self.origin = (model.mesh_xb[0], model.mesh_xb[2], model.mesh_xb[4])
        self._solid = None

    def center(self, i, j, k):
        ox, oy, oz = self.origin
        return (ox + (i + 0.5) * self.dx,
                oy + (j + 0.5) * self.dy,
                oz + (k + 0.5) * self.dz)

    def _build(self):
        nx, ny, nz = self.nx, self.ny, self.nz
        sol = bytearray(nx * ny * nz)

        def stamp(x0, x1, y0, y1, z0, z1, val):
            """把中心落在 [x0,x1]×[y0,y1]×[z0,z1] 内的格子置为 val（1=实心, 0=气）。"""
            i0 = max(0, int(math.floor((x0 - self.origin[0]) / self.dx)) - 1)
            i1 = min(nx, int(math.ceil((x1 - self.origin[0]) / self.dx)) + 1)
            j0 = max(0, int(math.floor((y0 - self.origin[1]) / self.dy)) - 1)
            j1 = min(ny, int(math.ceil((y1 - self.origin[1]) / self.dy)) + 1)
            k0 = max(0, int(math.floor((z0 - self.origin[2]) / self.dz)) - 1)
            k1 = min(nz, int(math.ceil((z1 - self.origin[2]) / self.dz)) + 1)
            for i in range(i0, i1):
                cx = self.origin[0] + (i + 0.5) * self.dx
                if not (x0 - 1e-9 <= cx <= x1 + 1e-9):
                    continue
                for j in range(j0, j1):
                    cy = self.origin[1] + (j + 0.5) * self.dy
                    if not (y0 - 1e-9 <= cy <= y1 + 1e-9):
                        continue
                    base = (i * ny + j) * nz
                    for k in range(k0, k1):
                        cz = self.origin[2] + (k + 0.5) * self.dz
                        if z0 - 1e-9 <= cz <= z1 + 1e-9:
                            sol[base + k] = val

        for b in self.m.obst:
            stamp(*b.xb, 1)
        for h in self.m.holes:          # &HOLE 把 OBST 挖成气相（0D/1D 反应器、窗洞等）
            stamp(*h, 0)
        self._solid = sol

    @property
    def solid(self):
        if self._solid is None:
            self._build()
        return self._solid

    def is_solid(self, i, j, k):
        if not (0 <= i < self.nx and 0 <= j < self.ny and 0 <= k < self.nz):
            return None          # 域外
        return bool(self.solid[(i * self.ny + j) * self.nz + k])

    def is_gas(self, i, j, k):
        return self.is_solid(i, j, k) is False

    def index_at(self, x, y, z):
        ox, oy, oz = self.origin
        i = int(math.floor((x - ox) / self.dx - 0.5 + 1e-9))
        j = int(math.floor((y - oy) / self.dy - 0.5 + 1e-9))
        k = int(math.floor((z - oz) / self.dz - 0.5 + 1e-9))
        return i, j, k

    def snap_error(self, v, axis):
        """面到最近网格平面的偏差，单位为该轴格子尺寸。"""
        ox = self.origin[axis]
        d = (self.dx, self.dy, self.dz)[axis]
        r = (v - ox) / d
        return abs(r - round(r))


# ──────────────────────────────────────────────────────────────────
# 检查
# ──────────────────────────────────────────────────────────────────

@dataclass
class Issue:
    code: str
    level: str      # 'ERROR' | 'WARN' | 'INFO'
    msg: str

    def __str__(self):
        mark = {"ERROR": "[X]", "WARN": "[!]", "INFO": "[i]"}[self.level]
        return f"{mark} {self.code:<18} {self.msg}"


def _inside_fraction(mesh, xb):
    """构件落在网格内的体积分数；退化轴按'平面是否在网格内'计。"""
    frac = 1.0
    for a in range(3):
        lo, hi = xb[2 * a], xb[2 * a + 1]
        mlo, mhi = mesh[2 * a], mesh[2 * a + 1]
        ext = hi - lo
        if ext < 1e-9:                                   # 零厚度面片
            frac *= 1.0 if mlo - 1e-6 <= lo <= mhi + 1e-6 else 0.0
        else:
            frac *= max(0.0, min(hi, mhi) - max(lo, mlo)) / ext
    return frac


def _inside_fraction_union(mesh_xbs, xb):
    """构件落在多网格**并集**内的体积分数（各轴独立）。退化轴按'面是否在并集内'计。

    ★ 不能用 max(单网格占比)：横跨多网格的构件（如全长地板）在每块网格里只占一部分，
      会误判成"被裁切"。必须把各网格在该轴的区间合并成并集后算覆盖长度。
    """
    frac = 1.0
    for a in range(3):
        lo, hi = xb[2 * a], xb[2 * a + 1]
        ext = hi - lo
        if ext < 1e-9:                                   # 零厚度面片
            inside = any(m[2 * a] - 1e-6 <= lo <= m[2 * a + 1] + 1e-6 for m in mesh_xbs)
            frac *= 1.0 if inside else 0.0
        else:
            ivs = sorted((max(lo, m[2 * a]), min(hi, m[2 * a + 1])) for m in mesh_xbs)
            covered = 0.0
            cur_lo = cur_hi = None
            for ilo, ihi in ivs:
                if ilo > ihi + 1e-12:
                    continue
                if cur_lo is None:
                    cur_lo, cur_hi = ilo, ihi
                elif ilo <= cur_hi + 1e-12:
                    cur_hi = max(cur_hi, ihi)
                else:
                    covered += cur_hi - cur_lo
                    cur_lo, cur_hi = ilo, ihi
            if cur_lo is not None:
                covered += cur_hi - cur_lo
            frac *= max(0.0, covered) / ext
    return frac


def _resolve_meshes(meshes):
    """把多网格平铺合并成单一逻辑网格（MPI 分解的常见形态）。

    返回 (ijk, xb, note)；无法合并（分辨率不一致 / 有间隙 / 有重叠）返回 None。
    均匀平铺网格的并集即原逻辑网格，合并后可用单网格逻辑精确校验。
    """
    if len(meshes) <= 1:
        return meshes[0][0], meshes[0][1], ""
    ijk0, xb0 = meshes[0]
    dx = (xb0[1] - xb0[0]) / ijk0[0]
    dy = (xb0[3] - xb0[2]) / ijk0[1]
    dz = (xb0[5] - xb0[4]) / ijk0[2]
    for ijk, xb in meshes:
        if abs((xb[1] - xb[0]) / ijk[0] - dx) > 1e-9 * max(1.0, abs(dx)):
            return None
        if abs((xb[3] - xb[2]) / ijk[1] - dy) > 1e-9 * max(1.0, abs(dy)):
            return None
        if abs((xb[5] - xb[4]) / ijk[2] - dz) > 1e-9 * max(1.0, abs(dz)):
            return None
    lo = [min(xb[2 * i] for _, xb in meshes) for i in range(3)]
    hi = [max(xb[2 * i + 1] for _, xb in meshes) for i in range(3)]
    nx = round((hi[0] - lo[0]) / dx)
    ny = round((hi[1] - lo[1]) / dy)
    nz = round((hi[2] - lo[2]) / dz)
    if sum(i * j * k for (i, j, k), _ in meshes) != nx * ny * nz:
        return None                      # 有间隙或重叠
    for (i, j, k), xb in meshes:
        if abs((xb[1] - xb[0]) / dx - i) > 1e-6:
            return None
        if abs((xb[3] - xb[2]) / dy - j) > 1e-6:
            return None
        if abs((xb[5] - xb[4]) / dz - k) > 1e-6:
            return None
    note = f"{len(meshes)} 个 &MESH 平铺合并为单网格 {nx}×{ny}×{nz}"
    return (nx, ny, nz), (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]), note


def check_model(m: FdsModel, resolution_floor: float = 4.0) -> list:
    out = []
    if not m.ijk or not m.mesh_xb:
        out.append(Issue("PARSE", "ERROR", "没解析到 &MESH IJK/XB，无法校验"))
        return out

    # 变换网格（TRNX/Y/Z_ID）无法按 Cartesian 近似校验，诚实跳过而非假报
    if m.transformed:
        out.append(Issue("MESH", "WARN",
                         "含变换网格(TRNX/Y/Z_ID)，本工具按 Cartesian 近似失效，跳过几何校验"))
        return out

    # 多网格：平铺合并成单一逻辑网格；非平铺则降级并明确告警（诚实限制）
    multi_note = ""
    if len(m.meshes) > 1:
        resolved = _resolve_meshes(m.meshes)
        if resolved is None:
            out.append(Issue("MESH", "WARN",
                             f"{len(m.meshes)} 个 &MESH 非均匀或含间隙/重叠，"
                             f"本工具按单一网格近似校验，结果可能不准"))
        else:
            m.ijk, m.mesh_xb, multi_note = resolved

    g = Grid(m)
    out.append(Issue("MESH", "INFO",
                     f"{m.chid or m.path.stem}: {m.cells:,} 格 "
                     f"({m.ijk[0]}×{m.ijk[1]}×{m.ijk[2]}, "
                     f"dx={m.dx:.4f}×{m.dy:.4f}×{m.dz:.4f} m)"))
    if multi_note:
        out.append(Issue("MESH", "INFO", multi_note))

    # ── E1/E2 出域与裁切 ──────────────────────────────────────
    # ★ 多网格/MULT：构件只在「不进任何一块网格」时才判域外；裁切按所有网格的
    #   **并集**体积占比算。不能用 max(单网格占比)：横跨多网格的构件（如全长
    #   地板）在每块网格里只占一部分，会误判成"被裁切"。
    mesh_xbs = [xb for _, xb in m.meshes] if m.meshes else [m.mesh_xb]
    outside, clipped = [], []
    for b in m.obst + m.vents:
        fr = _inside_fraction_union(mesh_xbs, b.xb)
        if fr < 1e-9:
            outside.append(b)
        elif fr < 0.999:
            clipped.append((b, fr))
    for b in outside:
        out.append(Issue("E1", "ERROR",
                         f"{b.kind} {b.name} 完全在域外 → FDS 静默丢弃。XB={b.xb}"))
    for b, fr in clipped:
        out.append(Issue("E2", "WARN",
                         f"{b.kind} {b.name} 被裁掉 {100*(1-fr):.0f}%。XB={b.xb} "
                         f"surf={b.surf}"))

    # ── G1 网格对齐 ──────────────────────────────────────────
    worst = []
    for b in m.obst + m.vents:
        e = max(max(g.snap_error(b.xb[2 * a], a), g.snap_error(b.xb[2 * a + 1], a))
                for a in range(3))
        if e > 1e-6:
            worst.append((e, b))
    if worst:
        worst.sort(key=lambda t: -t[0])
        out.append(Issue("G1", "WARN" if worst[0][0] > 0.01 else "INFO",
                         f"{len(worst)} 个构件的面不在网格平面上，"
                         f"最大偏移 {worst[0][0]:.3f} 格（FDS 会吸附到最近格子面）"))
        for e, b in worst[:6]:
            out.append(Issue("G1", "INFO",
                             f"  · {b.kind} {b.name} 偏移 {e:.3f} 格  XB={b.xb}"))

    # ── S1 吸附报告：declared → snapped（离散化程度） ─────────
    if worst:
        dd = (g.dx, g.dy, g.dz)
        out.append(Issue("S1", "INFO",
                         f"{len(worst)} 个面会被 FDS 静默吸附到最近格子面"))
        for e, b in worst[:6]:
            snapped = tuple(round(g.origin[a] + _nint((b.xb[2 * a + i] - g.origin[a]) / dd[a]) * dd[a], 4)
                            for a in range(3) for i in (0, 1))
            out.append(Issue("S1", "INFO",
                             f"  · {b.kind} {b.name}: {b.xb} → {snapped} "
                             f"(Δ={e:.3f} 格)"))

    # ── V1/V2 开口有没有气相通道 ─────────────────────────────
    # ★ 关键：VENT 落在网格边界上时**只有域内一侧存在格子**，另一侧是域外，
    #   不能当成"实心"记 0% —— 这是判据里最容易出假阳性的地方。
    d3 = (g.dx, g.dy, g.dz)
    n3 = (g.nx, g.ny, g.nz)
    for v in m.vents:
        a = v.is_plane()
        if not isinstance(a, int):
            out.append(Issue("V1", "WARN",
                             f"VENT {v.brief()} 不是单平面（退化轴 {a}），跳过通道检查"))
            continue
        pos = v.xb[2 * a]
        others = [x for x in range(3) if x != a]

        # 面片被吸附到哪条格子面（Fortran NINT 约定，半格向远离零）
        r = (pos - g.origin[a]) / d3[a]
        k = _nint(r)
        off = abs(r - k)

        # ★ 只统计**格心落在开口范围内**的格子。格心在开口外的格子不属于这个开口，
        #   把它算进来会得到假的"堵格"——门 z∈[0,2.2]、格心 2.3 的那一行在过梁里，
        #   但它压根不是门的一部分（曾把门的通流率误报成 92%）。
        rng = []
        for x in others:
            lo, hi = v.xb[2 * x], v.xb[2 * x + 1]
            i0 = int(math.ceil((lo - g.origin[x]) / d3[x] - 0.5 - 1e-9))
            i1 = int(math.floor((hi - g.origin[x]) / d3[x] - 0.5 + 1e-9))
            rng.append((max(0, i0), min(n3[x] - 1, i1)))
        cells = [(i, j) for i in range(rng[0][0], rng[0][1] + 1)
                 for j in range(rng[1][0], rng[1][1] + 1)]
        n_cell = len(cells)

        tag = f"VENT {v.brief()} 面 x[a={a}]={pos:.3f}"
        snap = f"，吸附偏移 {off:.2f} 格" if off > 1e-6 else ""

        if n_cell == 0:
            out.append(Issue("V1", "WARN",
                             f"{tag} 落格后没有任何格心落在开口内 → 该开口不产生流动"
                             f"。XB={v.xb}{snap}"))
            continue

        sides = {}
        for side, idx in (("下侧", k - 1), ("上侧", k)):
            if not (0 <= idx < n3[a]):
                sides[side] = (None, idx)          # 域外，不计
                continue
            n_gas = 0
            for i, j in cells:
                co = [0, 0, 0]
                co[a] = idx
                co[others[0]] = i
                co[others[1]] = j
                if g.is_gas(*co):
                    n_gas += 1
            sides[side] = (n_gas, idx)

        valid, dead = [], []
        for side, (n_gas, idx) in sides.items():
            (valid if n_gas is not None else dead).append(
                (side, n_gas, idx) if n_gas is not None else side)
        cov = " / ".join(f"{s} {n/n_cell:.0%}" for s, n, _ in valid)
        gone = f"（{'、'.join(dead)} 为域外，不计）" if dead else ""

        if valid and all(n == 0 for _, n, _ in valid):
            out.append(Issue("V1", "ERROR",
                             f"{tag} 在域内{'两侧' if len(valid) > 1 else '一侧'}全是实心格子"
                             f" → 该开口完全失效（FDS 静默忽略）。XB={v.xb}{snap}"))
        else:
            lv = "WARN" if (off > 0.01 and v.surf.upper() == "OPEN") else "INFO"
            out.append(Issue("V1", lv,
                             f"{tag} 气相覆盖 {cov}{gone}{snap}"))

        # ── V2 生效面积 vs 标称面积（部分堵口，仅限落格开口） ──
        # ★ 只对"面内尺寸已落格"的开口算生效面积。半格开口被 FDS 吸附是**平移**不是
        #   堵口——用声明尺寸数格心会把"偏移"误算成"部分被堵"，那是 G1/S1 的职责。
        # ★ 判据 = 某个域内侧 0 < 气相占比 < 95%（部分被堵）。全 0% 是 V1（全死），
        #   全 100% 是正常；火源 VENT 贴地板的那侧 0% 属"结构面"不报警（0 不进区间）。
        if valid:
            inplane_off = max(
                max(g.snap_error(v.xb[2 * x], x), g.snap_error(v.xb[2 * x + 1], x))
                for x in others)
            if inplane_off <= 0.01:
                ext = [v.xb[2 * o + 1] - v.xb[2 * o] for o in others]
                nom = ext[0] * ext[1]
                for side, n, _ in valid:
                    frac = n / n_cell
                    if 0.0 < frac < 0.95:
                        out.append(Issue("V2", "WARN",
                                         f"{tag} 标称面积 {nom:.3f} m²，{side} 生效仅 {frac:.0%}"
                                         f"（部分失效）"))
                        break

    # ── O1 重合/掩埋（盒子级） ───────────────────────────────
    obs = m.obst
    for a_i in range(len(obs)):
        for a_j in range(a_i + 1, len(obs)):
            A, B = obs[a_i], obs[a_j]
            if all(abs(A.xb[t] - B.xb[t]) < 1e-9 for t in range(6)):
                out.append(Issue("O1", "WARN",
                                 f"OBST {A.name} 与 {B.name} 完全重合（重复实体）"))
                continue
            if all(B.xb[t] - 1e-9 <= A.xb[t] and A.xb[t] <= B.xb[t] + 1e-9
                   for t in range(6)):
                out.append(Issue("O1", "WARN",
                                 f"OBST {A.name} 完全埋在 {B.name} 内部（几何无效）"))

    # ── O2 构件被其它构件吃掉的比例（格子级） ────────────────
    # ★ O1 只判"盒子完全包含"，抓不到这一类：C3 的电缆桥架 OBST 与侧墙在 y 上
    #   逐字相同、只在 z 上不同 —— 两者是**重叠**不是包含，O1 放行；
    #   但桥架占的每一格都已经被侧墙占了，桥架对整个流场零贡献。
    #   判据：某个 OBST 覆盖的格子里，有多少格在**别的** OBST 里也是实心。
    #   ≥99% → 该构件是惰性的（既不挡流也不传热），设计变量若挂在它身上就是死的。
    boxes_occl = []
    for idx_a, A in enumerate(obs):
        i0 = max(0, int(math.floor((A.xb[0] - g.origin[0]) / g.dx)) - 1)
        i1 = min(g.nx, int(math.ceil((A.xb[1] - g.origin[0]) / g.dx)) + 1)
        j0 = max(0, int(math.floor((A.xb[2] - g.origin[1]) / g.dy)) - 1)
        j1 = min(g.ny, int(math.ceil((A.xb[3] - g.origin[1]) / g.dy)) + 1)
        k0 = max(0, int(math.floor((A.xb[4] - g.origin[2]) / g.dz)) - 1)
        k1 = min(g.nz, int(math.ceil((A.xb[5] - g.origin[2]) / g.dz)) + 1)
        mine = tot = shared = 0
        for i in range(i0, i1):
            cx = g.origin[0] + (i + 0.5) * g.dx
            if not (A.xb[0] - 1e-9 <= cx <= A.xb[1] + 1e-9):
                continue
            for j in range(j0, j1):
                cy = g.origin[1] + (j + 0.5) * g.dy
                if not (A.xb[2] - 1e-9 <= cy <= A.xb[3] + 1e-9):
                    continue
                for k in range(k0, k1):
                    cz = g.origin[2] + (k + 0.5) * g.dz
                    if not (A.xb[4] - 1e-9 <= cz <= A.xb[5] + 1e-9):
                        continue
                    tot += 1
                    mine += 1
                    co = (cx, cy, cz)
                    for idx_b, B in enumerate(obs):
                        if idx_b == idx_a:
                            continue
                        # 第 t 轴要跟 xb[2t] / xb[2t+1] 比，不是 xb[t]
                        if all(B.xb[2 * t] - 1e-9 <= co[t] <= B.xb[2 * t + 1] + 1e-9
                               for t in range(3)):
                            shared += 1
                            break
        if mine:
            boxes_occl.append((shared / mine, A, mine))
    boxes_occl.sort(key=lambda t: -t[0])
    dead = [(r, b, n) for r, b, n in boxes_occl if r >= 0.99]
    # 被完全埋住 → 真缺陷（ERROR）当且仅当它是「燃烧/释气面」（火焰烧不到）
    # 或「真实材质」（带 MATL_ID/THICKNESS 的实体对象被埋，如纸3电缆埋柜体）；
    # 纯 INERT/默认面或显式 PERMIT_HOLE=F 的嵌套实体多为多层材料（钢芯、惰性
    # 隔离层），通常有意 → WARN。
    def _o2_is_defect(b):
        if b.permit_hole:
            return False
        return b.surf in m.burning_surfs or b.surf in m.material_surfs

    burn = [(r, b, n) for r, b, n in dead if _o2_is_defect(b)]
    inert = [(r, b, n) for r, b, n in dead if not _o2_is_defect(b)]
    if burn:
        out.append(Issue("O2", "ERROR",
                         f"{len(burn)} 个**燃烧/释气面** OBST 所覆盖的格子全部已被其它构件占用 → "
                         f"火焰无法触及，整个流场零贡献（设计变量挂在这些构件上等于死变量）"))
        for r, b, n in burn[:12]:
            out.append(Issue("O2", "ERROR",
                             f"  · {b.name}(surf={b.surf}) {n} 格 100% 被占  XB={b.xb}"))
        if len(burn) > 12:
            out.append(Issue("O2", "ERROR", f"  … 另有 {len(burn)-12} 个"))
    if inert:
        out.append(Issue("O2", "WARN",
                         f"{len(inert)} 个惰性 OBST 被完全包在其它实体内（多层材料/嵌套实体，通常有意）"))
        for r, b, n in inert[:6]:
            out.append(Issue("O2", "WARN",
                             f"  · {b.name}(surf={b.surf}) {n} 格 100% 被占  XB={b.xb}"))
    # 只报 [0.9, 0.99) 这一档：≥0.99 已是 ERROR，0.5~0.9 大多是**合法重叠**
    # （例如侧墙的格子被贴在墙上的托盘共用），报出来只会淹没真信号。
    part = [(r, b, n) for r, b, n in boxes_occl if 0.9 <= r < 0.99]
    if part:
        out.append(Issue("O2", "WARN",
                         f"{len(part)} 个 OBST 半数以上格子被其它构件占用"))
        for r, b, n in part[:6]:
            out.append(Issue("O2", "WARN",
                             f"  · {b.name} 被占 {r:.0%}（{n} 格中）  XB={b.xb}"))

    # ── D1 传感器落在实心里 ──────────────────────────────────
    # ★ 边界点（如 z=0、y=0）会被 index_at 的最近格心公式映到 -1，但那是**域内边界**
    #   不是域外。先按物理坐标判"真在域外"，再钳位到最近格子判实心。
    for did, q, xyz in m.devc:
        if not any(_point_in_xb(xyz, mx) for mx in mesh_xbs):
            out.append(Issue("D1", "ERROR",
                             f"DEVC '{did}' ({q}) XYZ={xyz} 落在网格外"))
            continue
        if (q or "").upper() in _WALL_Q | _NONSPATIAL_Q:
            continue                       # 墙面/计算量落在实心是合法用法
        i = max(0, min(g.nx - 1, g.index_at(*xyz)[0]))
        j = max(0, min(g.ny - 1, g.index_at(*xyz)[1]))
        k = max(0, min(g.nz - 1, g.index_at(*xyz)[2]))
        if g.is_solid(i, j, k) is True:
            out.append(Issue("D1", "ERROR",
                             f"DEVC '{did}' ({q}) XYZ={xyz} 落在实心格子里，读数无意义"))

    # ── B1 网格边界退化审计 ──────────────────────────────────
    # 判据：网格边界上若出现**气相格子**，该格子所在的外表面就不再是墙体，
    # 而退化成 FDS 默认的 INERT 面。被 VENT（门/排烟口）覆盖的属于设计开口，
    # 剩下的才是"本该是混凝土、实际消失了"的假墙。
    d3 = (g.dx, g.dy, g.dz)
    n3 = (g.nx, g.ny, g.nz)
    covers = []      # (面名, 气相格数, 被 VENT 覆盖的格数)

    def _face_plane(axis, which):
        """which=0 → 下边界平面；which=1 → 上边界平面。"""
        if which == 0:
            return g.origin[axis]
        return g.origin[axis] + n3[axis] * d3[axis]

    def _vent_covers(v, axis, which, i, j, k):
        """该 VENT 是否盖住这个边界格子（平面同轴、同侧、横向范围内）。"""
        if v.is_plane() != axis:
            return False
        pos = v.xb[2 * axis]
        if abs(pos - _face_plane(axis, which)) > 1e-6:
            return False
        p = list(g.center(i, j, k))
        p[axis] = pos
        return all(v.xb[2 * a] - 1e-9 <= p[a] <= v.xb[2 * a + 1] + 1e-9
                   for a in range(3))

    faces = [("x-min", 0, 0), ("x-max", 0, 1),
             ("y-min", 1, 0), ("y-max", 1, 1),
             ("z-min", 2, 0), ("z-max", 2, 1)]
    for fname, axis, which in faces:
        idx = 0 if which == 0 else n3[axis] - 1
        gas = cov = 0
        for i in range(n3[0]):
            for j in range(n3[1]):
                for k in range(n3[2]):
                    if (i, j, k)[axis] != idx:
                        continue
                    if not g.is_gas(i, j, k):
                        continue
                    gas += 1
                    if any(_vent_covers(v, axis, which, i, j, k) for v in m.vents):
                        cov += 1
        covers.append((fname, gas, cov))

    tot_gas = sum(c[1] for c in covers)
    if tot_gas:
        bare = tot_gas - sum(c[2] for c in covers)
        detail = "  ".join(f"{n}={gg}(开口{cv},裸{gg-cv})"
                          for n, gg, cv in covers if gg)
        lv = "WARN" if bare > 0 else "INFO"
        out.append(Issue("B1", lv,
                         f"网格边界气相格 {tot_gas} 个：{detail}"))
        if bare > 0:
            out.append(Issue("B1", lv,
                             f"     其中 {bare} 格（{bare/tot_gas:.0%}）无 VENT 覆盖 → "
                             f"本该是混凝土墙，实际退化成 FDS 默认 INERT 面"))

    # ── R1 网格分辨率 ────────────────────────────────────────
    if m.hrrpua_sum_area > 0:
        rho, cp, T, grav = 1.204, 1.005, 293.15, 9.81
        q = m.hrrpua_sum_area
        dstar = (q / (rho * cp * T * math.sqrt(grav))) ** 0.4
        ratio = dstar / min(m.dx, m.dy, m.dz)
        lv = "INFO" if ratio >= resolution_floor else "WARN"
        out.append(Issue("R1", lv,
                         f"火源 ΣHRRPUA·A ≈ {q:.0f} kW → D*={dstar:.3f} m, "
                         f"D*/dx={ratio:.2f}"
                         + ("" if ratio >= resolution_floor else
                            f"  < {resolution_floor:g}（低于常用准则）")))

    return out


def check_file(path: Path, quiet_ok: bool = False):
    issues = check_model(parse_fds(path))
    if quiet_ok and not any(i.level == "ERROR" for i in issues):
        return issues
    print(f"\n{'='*78}\n{path}\n{'='*78}")
    for i in issues:
        print("  " + str(i).replace("\n", "\n  "))
    return issues


def report(paths):
    n_err = n_warn = n_files = 0
    bad = []
    for p in paths:
        iss = check_model(parse_fds(p))
        n_files += 1
        e = sum(1 for i in iss if i.level == "ERROR")
        w = sum(1 for i in iss if i.level == "WARN")
        n_err += e
        n_warn += w
        if e:
            bad.append((p, e, w))
    print(f"\n共校验 {n_files} 个 .fds：ERROR {n_err} 条 / WARN {n_warn} 条")
    if bad:
        print("\n有 ERROR 的文件：")
        for p, e, w in bad:
            print(f"  {e:>3} err  {w:>3} warn   {p}")
    return n_err, n_warn


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if not args:
        print(__doc__)
        return 1
    targets = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            targets += sorted(p.rglob("*.fds"))
        elif p.exists():
            targets.append(p)
        else:
            print(f"跳过（不存在）：{a}")
    if len(targets) == 1:
        check_file(targets[0])
    else:
        report(targets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
