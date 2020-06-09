import calkulate as calk


def test_authors():
    assert isinstance(calk.__authorlist__, list)
    assert isinstance(calk.__author__, str)


def test_version():
    assert isinstance(calk.__version__, str)
