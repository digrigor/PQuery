import os
sqldb = 'sqlite:///app/resources/varDB'


class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'

SQLALCHEMY_DATABASE_URI = sqldb
SQLALCHEMY_TRACK_MODIFICATIONS = False
