from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32

with open('../../collector.recordio', 'rb') as f:
    buf = f.read()
    n = 0
    while n < len(buf):
        msg_len, new_pos = _DecodeVarint32(buf, n)
        n = new_pos
        msg_buf = buf[n:n+msg_len]
        n += msg_len
        read_metric = metric_pb2.Metric()
        read_metric.ParseFromString(msg_buf)
        # do something with read_metric

# $ protoc -I=$SRC_DIR --python_out=$DST_DIR $SRC_DIR/output.proto
# -> generates output.py
