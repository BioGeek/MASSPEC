from masspec.example import pxd071026


def test_packaged_run_loads():
    run = pxd071026()
    assert len(run.masses) == 5440
    assert len(run.myoglobin) == 154 and "I" not in run.myoglobin
    assert run.pepsin.startswith("EPLENYLDTEYFGT")
