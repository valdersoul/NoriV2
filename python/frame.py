import numpy as np
import math


def sinTheta2(v):
    return 1.0 - v[2] * v[2]


def sinTheta(v):
    return safe_sqrt(sinTheta2(v))


def cosTheta(v):
    return v[2]


def tanTheta(v):
    return sinTheta(v) / cosTheta(v)

def safe_sqrt(value):
    return math.sqrt(max(0.0, value))
