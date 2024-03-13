#!/bin/bash

flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
flake8 . --count --exit-zero --max-complexity=13 --max-line-length=127 --per-file-ignores="__init__.py:F401" --statistics
