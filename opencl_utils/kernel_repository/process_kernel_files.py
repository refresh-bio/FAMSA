import sys
import os


className = 'KernelRepository'
dir = '.'

classDecl = open(className + '.h', 'w')
classDecl.write(
	'#pragma once\n\n'
	'class ' + className + '\n'
	'{\n'
	'public:\n'
	'static char * kernels[];\n')

classDef = open(className + '.cpp', 'w')
classDef.write('#include \"' + className + '.h\"\n\n')

counter = 0
fileList = []

for fileName in os.listdir(dir): 
	if os.path.isfile(os.path.join(dir,fileName)):
		splitted = os.path.splitext(fileName)
		base = splitted[0]
		extension = splitted[1]
		
		if extension == '.cl':
			fileList.append(fileName)
					
			classDecl.write('static const int ' + base + ';\n')
			classDef.write('const int ' + className + '::' + base + ' = ' + str(counter) + ';\n')
			counter = counter + 1
			
classDef.write('char * ' + className + '::kernels[] = {\n')
		
counter = 0
for fileName in fileList:	
	file = open(fileName, 'r')	
	for line in file:
		if line[:3] == '///':
			continue;
	
		line = line.strip()
		line = line.replace('"', '\\"');
		line = line.replace('\\n', '\\\\n');
		
		classDef.write('\"' + line + '\\n\"\n')
		
	classDef.write(',\n')
	counter = counter + 1
				
classDecl.write('};')			
classDef.write('""\n};')