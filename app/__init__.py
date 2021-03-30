from flask import Flask, jsonify, render_template, request
from config import Config
import pandas as pd
import sqlite3
from datatables import ColumnDT, DataTables
from flask_sqlalchemy import SQLAlchemy
from app.handlers import *
from app.dependencies import *

app = Flask(__name__)
app.secret_key = 'A0Zr98j/3yX R~XHH!jmN]LWX/,?RT'
app.config.from_object(Config)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///resources/poster_072019.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False


#comm = 'dt_columns = [\n' + ',\n'.join(['ColumnDT(Var.' + x + ')' for x in db_cols]) + '\n]'
#exec(comm)

#print(filter_dict)
from app import routes
