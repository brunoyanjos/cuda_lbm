def find_idx(list, ref):
    for idx, value in enumerate(list):
        if value > ref:
            return [idx - 1, idx]
        elif value == ref:
            return [idx]
    
    return []

def tke_interpolation(t1, tke1, t2, tke2, t3):
    return tke1 + (tke2 - tke1) * (t3 - t1)/(t2 - t1)