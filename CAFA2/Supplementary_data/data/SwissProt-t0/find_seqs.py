from pymongo import MongoClient
mongoCli = MongoClient("localhost")
mongoDB = mongoCli.eggnog4_1
db_members = mongoDB.members

db_members
