import sys

if len(sys.argv) != 2 :
  sys.stderr.write('Usage: {} <nucleotide file>\n'.format(sys.argv[0]))
  sys.exit()

# assumes the file only contains dna and newlines
newbytearray=bytearray(b'',encoding='utf-8')
dict={'A':0b1000,'C':0b0100,'G':0b0010,'T':0b0001,'\n':0b1010}
with open(sys.argv[1]) as file:
    while True:
        char=file.read(1)
        if not char:
            file.close()
            break
        newbytearray.append(dict[char])
outfile = open(sys.argv[1] + '.bin', 'wb')
outfile.write(newbytearray)
outfile.close()

#Converts the binary file to unicode and prints the result sequence.
testBin = open('fileA.txt.bin','rb')
sequence=''
for line in testBin:
    line = line.replace(chr(0b1000),'1000')
    line = line.replace(chr(0b0100),'0100')
    line = line.replace(chr(0b0010),'0010')
    line = line.replace(chr(0b0001),'0001')
    line = line.replace(chr(0b1010),'\n')
    sequence += line
#outputVerify = open('outputVerify.txt','wb')
#outputVerify.write(sequence)
#outputVerify.close()
print sequence
testBin.close()

#Shows the data of the binary file. Note that byte 6 is the newline character 0b1010.
testBin = open('fileA.txt.bin','rb')
list = ''
i=0
while True:
    b = testBin.read(1)
    i += 1
    if not b:
    break #due to eof
    list += b
    print 'byte: ' + str(i) + ' is '+ '{0:04b}'.format(ord(b)) +' and has decimal representation: ' + str(ord(b))
testBin.close()
