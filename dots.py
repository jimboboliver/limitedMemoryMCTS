import os
from subprocess import call
from sys import argv

if len(argv) > 1 :
    call(["dot", argv[1], "-Tjpg", "-o", argv[1].split('.')[0] + '.jpg'])
else:
    for item in os.listdir('graph'):
        if 'graph' in item and '.dot' in item:
            call(["dot", 'graph/' + item, "-Tjpg", "-o", 'graph/' + item.split('.')[0] + '.jpg'])
