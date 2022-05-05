# -*- coding: utf-8 -*-
"""
Created on Tue May  3 15:21:20 2022

@author: tdeguire
"""


cf = open('tape41', 'r') #tape[i]
test = []

for line in cf:
    test.append(line.split())
    

cf.close()