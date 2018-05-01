import os
from subprocess import call

for item in os.listdir('.'):
    if 'graph' in item and '.dot' in item:
        call(["dot", item, "-Tjpg", "-o", item.split('.')[0] + '.jpg'])
