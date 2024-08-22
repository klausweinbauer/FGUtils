#!/bin/bash

flake8 . --count --exclude venv/ --select=E9,F63,F7,F82 --show-source --statistics
flake8 . --count --exclude venv/ --exit-zero --max-complexity=13 --max-line-length=127 --per-file-ignores="__init__.py:F401,fgutils/parse.py:E203" --statistics
