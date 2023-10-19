# Some datasets have the non consecutive edges in the wrong place. This script fixes the dataset

filename = "parking-garage.g2o"

f = open(filename, 'r')
lines = f.readlines()

# stores non consecutive edges
d = {}
# stores consecutive edges
contig = {}
# vertices
vert = ''

for line in lines:
    print(line)
    split_line = line.split(' ')
    if 'VERT' in split_line[0]:
        vert += line
        continue
    key = max(int(split_line[1]), int(split_line[2]))
    if abs(int(split_line[1]) - int(split_line[2])) == 1:
        contig[key] = line
        continue
    if key not in d:
        d[key] = ''
    d[key] += line

f.close()

f = open("parking-garage-fixed.g2o", 'w')

f.write(vert)
for key in contig:
    f.write(contig[key])
    if key in d:
        f.write(d[key])

f.close()
    

    
