from timeit import timeit

from pkdb_analysis import PKData

from tests import TEST_PKDATA


def bla():
    PKData.from_archive(TEST_PKDATA)
