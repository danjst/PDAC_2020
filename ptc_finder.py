f = open("plus_131.fa")

dna_dict = {}
sequence = ""
name = ""

for line in f:
    if line[0] == ">":
        if name != "":
            dna_dict[name] = sequence
        sequence = ""
        name = line.split(" ")[0][1:]
    else:
        sequence += line.rstrip("\n")
dna_dict[name] = sequence

f = open("output.csv", "w")
output_arr = {}
for x in dna_dict:
    sequence = dna_dict[x]
    i = 0
    found = False
    dna_found = ""
    while len(sequence) > 3 and not found:
        triple = sequence[0:3]
        if triple == "TAA":
            dna_found = "TAA"
            found = True
        elif triple == "TAG":
            dna_found = "TAG"
            found = True
        elif triple == "TGA":
            dna_found = "TGA"
            found = True
        sequence = sequence[3:]
        i = i + 1
    if found:
        output_arr[x] = x + "," + str(dna_found) + ":" + str(3 * (i - 1) + 1)
    else:
        output_arr[x] = x + ",not found"

for x in dna_dict:
    sequence = dna_dict[x][1:]
    i = 0
    found = False
    dna_found = ""
    while len(sequence) > 3 and not found:
        triple = sequence[0:3]
        if triple == "TAA":
            dna_found = "TAA"
            found = True
        elif triple == "TAG":
            dna_found = "TAG"
            found = True
        elif triple == "TGA":
            dna_found = "TGA"
            found = True
        sequence = sequence[3:]
        i = i + 1
    if found:
        output_arr[x] += "," + str(dna_found) + ":" + str(3 * (i - 1) + 2)
    else:
        output_arr[x] += ",not found"

for x in dna_dict:
    sequence = dna_dict[x][2:]
    i = 0
    found = False
    dna_found = ""
    while len(sequence) > 3 and not found:
        triple = sequence[0:3]
        if triple == "TAA":
            dna_found = "TAA"
            found = True
        elif triple == "TAG":
            dna_found = "TAG"
            found = True
        elif triple == "TGA":
            dna_found = "TGA"
            found = True
        sequence = sequence[3:]
        i = i + 1
    if found:
        output_arr[x] += "," + str(dna_found) + ":" +  str(3 * (i - 1) + 3) + "\n"
    else:
        output_arr[x] += ",not found\n"

f.write("Name,Frame1,Frame2,Frame3\n")
for dna in output_arr:
    f.write(output_arr[dna])
