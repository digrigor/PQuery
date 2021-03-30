#!/bin/bash

source venv_py3/bin/activate
export FLASK_ENV=development
export FLASK_APP=PQuery.py
flask run
