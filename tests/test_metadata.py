import calkulate as calk


def test_authors():
    assert isinstance(calk.__author__, str)


def test_version():
    assert isinstance(calk.__version__, str)
    calk.say_hello()


test_authors()
test_version()
