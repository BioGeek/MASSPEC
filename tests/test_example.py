from masspec.example import pxd071026, pxd071026_instanovo


def test_packaged_run_loads():
    run = pxd071026()
    assert len(run.masses) == 5440
    assert len(run.myoglobin) == 154 and "I" not in run.myoglobin
    assert run.pepsin.startswith("EPLENYLDTEYFGT")


def test_packaged_predictions_load():
    predictions = pxd071026_instanovo()
    assert len(predictions) == 3843
    assert all(0.0 <= row.confidence <= 1.0 for row in predictions)
