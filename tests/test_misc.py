# %%
import calkulate as calk


def test_hello():
    # Say hello
    calk.hello()
    calk.say_hello()


def test_metadata():
    """Are the metadata strings there?"""
    assert isinstance(calk.__version__, str)
    assert isinstance(calk.__author__, str)


# test_hello()
# test_metadata()
