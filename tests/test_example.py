from masspec.example import (
    acid_hydrolysates,
    pxd071026,
    pxd071026_instanovo,
)


def test_packaged_run_loads():
    run = pxd071026()
    assert len(run.masses) == 5440
    assert len(run.myoglobin) == 154 and "I" not in run.myoglobin
    assert run.pepsin.startswith("EPLENYLDTEYFGT")


def test_packaged_predictions_load():
    predictions = pxd071026_instanovo()
    assert len(predictions) == 3843
    assert all(0.0 <= row.confidence <= 1.0 for row in predictions)


def test_packaged_hydrolysates_load():
    runs = acid_hydrolysates()
    assert [run.accession for run in runs] == ["PXD003804", "PXD030203"]
    assert [len(run.reads) for run in runs] == [5449, 16462]
    for run in runs:
        assert len(run.myoglobin) == 153
        assert run.myoglobin.startswith("GLSDGEWQQVLNVWGKVEADLAGHGQ")
        assert all(read.alc >= 50 for read in run.reads)
