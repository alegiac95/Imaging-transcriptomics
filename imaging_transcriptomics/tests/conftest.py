import pytest
import  numpy as np
from pathlib import Path


@pytest.fixture(scope="session")
def tdata_dir():
    yield Path(__file__).parent / 'data'

