#!/bin/bash

PYTHONPATH=../../python OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python train.py
#PYTHONPATH=../../python python train.py
