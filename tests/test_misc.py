import calkulate as calk

# Say hello
calk.hello()


def test_metadata():
    """Are the metadata strings there?"""
    assert isinstance(calk.__version__, str)
    assert isinstance(calk.__author__, str)


# test_metadata()
