import pytest
from rvic.core.log import StreamToFile


def test_setup_stream():
    stream = StreamToFile()
    stream.write("buf")
    stream.flush()


def test_stream_raises_with_bad_log_level():
    with pytest.raises(TypeError):
        stream = StreamToFile(log_level="junk")
        stream.write("buf")
        stream.flush()


# Cannot test logger in interactive session or using pytest
