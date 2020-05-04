import numpy as np
from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32
from pp_mix.protos.py.state_pb2 import EigenMatrix, EigenVector


def gen_even_slices(n, n_packs, n_samples=None):
    start = 0
    if n_packs < 1:
        raise ValueError("gen_even_slices got n_packs=%s, must be >=1"
                         % n_packs)
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            if n_samples is not None:
                end = min(n_samples, end)
            yield np.arange(start, end)
            start = end


def loadChains(filename, msgType):
    out = []
    with open(filename, "rb") as fp:
        buf = fp.read()

    n = 0
    while n < len(buf):
        msg_len, new_pos = _DecodeVarint32(buf, n)
        n = new_pos
        msg_buf = buf[n:n+msg_len]
        try:
            msg = msgType()
            msg.ParseFromString(msg_buf)
            out.append(msg)
            n += msg_len
        except Exception as e:
            break

    return out


def writeChains(chains, filename):
    with open(filename, "wb") as fp:
        for c in chains:
            try:
                msgStr = c.SerializeToString()
                delimiter = _VarintBytes(len(msgStr))
                fp.write(delimiter + msgStr)
            except Exception as e:
                print(e)
                break   


def to_numpy(obj):
    if isinstance(obj, EigenMatrix):
        out = np.array(obj.data).reshape(obj.rows, obj.cols, order="F")
    elif isinstance(obj, EigenVector):
        out = np.array(obj.data)
    else:
        raise ValueError("incorrect object type")

    return out


def to_proto(array):
    if array.ndim == 1:
        out = EigenVector()
        out.size = len(array)
        out.data.extend(array.tolist())
    else:
        out = EigenMatrix()
        out.rows = array.shape[0]
        out.cols = array.shape[1]
        out.data.extend(array.reshape(1, -1, order='F').tolist()[0])
    return out