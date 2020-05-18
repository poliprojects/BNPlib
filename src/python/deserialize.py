from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32
import output_pb2
from collections import deque

with open('../../collector.recordio', 'rb') as f:
    buf = f.read()
    n = 0
    d=deque()
    while n < len(buf):
        msg_len, new_pos = _DecodeVarint32(buf, n)
        n = new_pos
        msg_buf = buf[n:n+msg_len]
        n += msg_len
        read_metric = output_pb2.State()
        read_metric.ParseFromString(msg_buf)
        d.append(read_metric)
    print(d[3])
      
