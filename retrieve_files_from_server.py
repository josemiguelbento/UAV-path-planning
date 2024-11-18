import paramiko
import os

client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect('', username='', password='')


sftp = client.open_sftp()

file_name = ''
folder = "results"

i = 0
for file_name in sftp.listdir('./UAV-path-planning/'+folder+'/'):
    if not os.path.isfile("./"+folder+"/"+file_name):
        sftp.get("./UAV-path-planning/"+folder+"/" + file_name, "./"+folder+"/"+file_name)
        print("./UAV-path-planning/"+folder+"/" + file_name, "./"+folder+"/"+file_name)
        i+=1
print(i," files were downloaded")
sftp.close()