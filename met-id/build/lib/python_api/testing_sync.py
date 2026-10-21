from sync_server import MetIDServer

server = MetIDServer(project_dir="../src-tauri", binary_name="")

server.start()

print(server.return_name("Patrik"))

data = [i for i in range(10)]
print(server.process_items(data))
try:
    print(server.identify())
    server.stop()
except:
    server.stop()

