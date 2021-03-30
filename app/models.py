from app import app
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy(app)
db.init_app(app)
db.Model.metadata.reflect(db.engine)

class Var(db.Model):
    __tablename__ = 'variant_info'
    __table_args__ = { 'extend_existing': True }
    variant_id = db.Column(db.Text, primary_key=True, index=True, unique=True)

class Alleles(db.Model):
    __tablename__ = 'variant_geno'
    __table_args__ = { 'extend_existing': True }
    variant_id = db.Column(db.Text, primary_key=True, index=True)

class FieldTypes(db.Model):
    __tablename__ = 'processed_field_types'
    __table_args__ = { 'extend_existing': True }
    index = db.Column(db.Text, primary_key=True)

class SelectValues(db.Model):
    __tablename__ = 'select_unique_values'
    __table_args__ = { 'extend_existing': True }
    index = db.Column(db.Text, primary_key=True)

class Samples(db.Model):
    __tablename__ = 'samples'
    __table_args__ = { 'extend_existing': True }
    name = db.Column(db.Text, primary_key=True)

class Temp(db.Model):
    __tablename__ = 't'
    __table_args__ = { 'extend_existing': True }
    vid = db.Column(db.Text, primary_key=True, index=True, unique=True)