"""Smoke tests for fds_audit."""

from fds_audit import check_model, parse_text


def _box(extra: str = "") -> str:
    return f"""
&HEAD CHID='test' /
&TIME T_END=10.0 /
&MESH IJK=20,20,20 XB=0.0,4.0,0.0,4.0,0.0,4.0 /
{extra}
"""


def test_parse_minimal():
    m = parse_text(_box())
    assert m.chid == "test"
    assert m.ijk == (20, 20, 20)


def test_clean_box_has_no_error():
    issues = check_model(parse_text(_box()))
    assert not any(i.level == "ERROR" for i in issues)


def test_e1_obst_outside_mesh():
    txt = _box("&OBST XB=4.5,5.0,1.0,3.0,0.0,2.0 SURF_ID='INERT' /\n")
    issues = check_model(parse_text(txt))
    assert any(i.code == "E1" and i.level == "ERROR" for i in issues)


def test_v1_dead_vent_between_solids():
    # A VENT sandwiched between two solid slabs has no gas on either side.
    txt = _box(
        "&OBST XB=0.0,4.0,0.0,4.0,0.0,0.2 SURF_ID='INERT' /\n"
        "&OBST XB=0.0,4.0,0.0,4.0,0.2,0.4 SURF_ID='INERT' /\n"
        "&VENT XB=1.0,3.0,1.0,3.0,0.2,0.2 SURF_ID='OPEN' /\n"
    )
    issues = check_model(parse_text(txt))
    assert any(i.code == "V1" and i.level == "ERROR" for i in issues)


def test_b1_mb_vent_not_bare():
    # Full-boundary &VENT MB='XMIN'...'OPEN' must not be flagged as bare wall (B1 WARN).
    txt = _box(
        "&VENT MB='XMIN', SURF_ID='OPEN' /\n"
        "&VENT MB='XMAX', SURF_ID='OPEN' /\n"
        "&VENT MB='YMIN', SURF_ID='OPEN' /\n"
        "&VENT MB='YMAX', SURF_ID='OPEN' /\n"
        "&VENT MB='ZMIN', SURF_ID='OPEN' /\n"
        "&VENT MB='ZMAX', SURF_ID='OPEN' /\n"
    )
    issues = check_model(parse_text(txt))
    b1 = [i for i in issues if i.code == "B1"]
    assert b1, "B1 应仍报告边界气相格统计"
    assert all(i.level != "WARN" for i in b1), "MB 风口已覆盖的边界不应有 B1 WARN"


def test_b1_mb_faces_parsed():
    m = parse_text(_box("&VENT MB='XMIN', SURF_ID='OPEN' /\n"))
    assert (0, 0) in m.mb_faces
    assert m.vents == []  # MB-only vent has no XB, so no Box is created
