#!/bin/sh
cp setup.py setup_tmp.py
cp setup_dev.py setup.py
pip install . -v
rm setup.py
mv setup_tmp.py setup.py