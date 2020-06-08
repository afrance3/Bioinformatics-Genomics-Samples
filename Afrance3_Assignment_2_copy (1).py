'''
Adam m France
Assignment 2

'''
aa_dict = {
'Met':['ATG'],
'Phe':['TTT', 'TTC'],
'Leu':['TTA', 'TTG', 'CTT','CTC', 'CTA', 'CTG'], 
'Cys':['TGT', 'TGC'], 
'Tyr':['TAC', 'TAT'], 
'Trp':['TGG'],
'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 
'His':['CAT', 'CAC'], 
'Gln':['CAA', 'CAG'],
'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
'Ile':['ATT', 'ATC', 'ATA'],
'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 
'Asn':['AAT', 'AAC'], 
'Lys':['AAA', 'AAG'],
'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 
'Val':['GTT', 'GTC', 'GTA', 'GTG'],
'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 
'Asp':['GAT', 'GAC'], 
'Glu':['GAA', 'GAG'],
'Gly':['GGT', 'GGC', 'GGA', 'GGG'],
'*':['TAA','TAG','TGA']
 }

#read in file
path = r"C:\Users\Adam_linux\Downloads\assignment2v2.fastq"
f = open(path,'r')
data = f.readlines() #file is now in list form

#remove newline characters?
cleanlist = []
for i in range( len(data)):
    cleanlist.append(data[i].rstrip('\n'))
'''
1. Search record header lines, and return records that match
the search criteria
'''
def function1(search,data):

    returnlist = []  #initialize list
    for i in range(0, len(data),4):  #loop through data to find terms
        if search in data[i]:   #if search is found return to list
            returnlist.append( data[i])
            returnlist.append( data[(i+1)])
            returnlist.append( data[(i+2)])
            returnlist.append( data[(i+3)])
    print('found '+str(len(returnlist)))
    print(returnlist)
    return returnlist
'''
2. Select a record, and generate a count of each amino acid coded for by the codons
of this sequence record.  Keep in mind that because these records are not necessarily
in the proper reading frame, so the user should be prompted to select a reading
frame (0, +1, +2).  Keep in mind that if you build your data structure of the FASTQ
file using a dictionary, that there is no order to the records in dictionaries.
How you plan to present the choice of records to the user during this step may
require some cleverness.  Presenting only a subset of the total 25 records might be a
good idea.
'''

def function2(cutoff,selection):   

    #sequence will be index 'select + 1'
    x = cleanlist[selection+1]
    
    seq = []
    #seq becomes a list of lists with DNA sequence
    if ask == 0:
        for i in range(0,len(x),3):
            seq.append(x[i:(i+3)])  #reading frame 0
    elif ask == 1:
        for i in range(0,len(x),3):
            seq.append(x[(i+1):(i+3+1)])#reading frame +1
    elif ask == 2:
        for i in range(0,len(x),3):
            seq.append(x[(i+2):(i+3+2)])#reading frame +2

    #iterate through seq and translate codons into amino acids
    d = dict()
    aalist = []
    for i in seq:
        for j in aa_dict:
            if i in aa_dict[j]:
                aalist.append(j)

    #count aalist into new d
    for i in aalist:
        d[i] = d.get(i,0) + 1

    print(d)
    return d 

'''
3. Allow the user to trim the nucleotide sequence record based on the quality score.
User should be prompted for a score cutoff, after which a new file
("assignment2_trimmed.fastq") should be generated of the FASTQ records with the sequence
and quality score lines trimmed based on that selection.
'''
def function3(score,cleanlist):

    scorelist = []#the quality scores from cleanlist will fo here
    for i in range(3, len(cleanlist),4):
        scorelist.append(cleanlist[i])

    qual_list=[]#same as scorelist, but in a list of lists format

    for i in range( len(scorelist)):
        qual_list.append(list(scorelist[i]))

    list3 = []#same as qual_list but in integer form
    for line in range( len(qual_list)): 
            list4=[]
            for i in qual_list[line]:
                x = ord(i)-64
                list4.append(x)
            list3.append(list4)

    cutoff=[]   #this will give us a list of the indicies where the avg score fell below the user's scpec
    for x in range( 0, len(list3)):
        count = 0
        for i in range( 0, len(list3[x])):
            avg = ( list3[x][i] + list3[x][i+1] + list3[x][i+2] ) / 3
            count += 1
            if avg >= score:
                continue
            else:
                cutoff.append(count)
                break
    #this bloc prepares list for export to the new file       
    outlist = cleanlist.copy() # a copy of the list form of our data
    count = 0
    for i in range( 1, len(cleanlist[i]), 4): # sequences
            x = ( cleanlist[i][: cutoff[count]])
            outlist[i] = x
            count += 1
    count = 0
    for i in range( 3, len(cleanlist[i]), 4): #quality scores
            x = ( cleanlist[i][: cutoff[count]])
            outlist[i] = x
            count += 1

        #outlist now ready to be wrritten to file
    out = open(r"C:\Users\Adam_linux\Downloads\assignment2_trimmed.fastq", 'w')
    for i in range( len(outlist)):
        out.write(outlist[i]+'\n')
    out.close()
    return outlist    

   

'''
These three functions should be called from a main method, that like the previous
assignment, collects user inputs and directions, then executes whatever the user
selects. The three functions should use the various user inputs and return a result.
The functions should not ask for user input themselves, they should only interpret
arguments and execute a task. 
'''


def main():
    #menu, not linear
    menu = int(input('please select from the following: \n1. Search record by keyword\n2. Amino acid count \n3. Sequence trim'))
    
    if menu == 1:
               search = input('Enter your search word: ')
               function1(search,cleanlist)
    
    elif menu == 2:
        for i in range(0,len(cleanlist),4):
            print(cleanlist[i])
        #have user select from the output
        select = int(input('please enter the number to left of the record that you want\n '))
        #select a record and generate a count of each amino acid codded for by the codons
        #of this sequence record.
        rframe = int(input('select a reading frame: \n 0\n 1\n 2\n'))
        function2(rframe,select)

    elif menu == 3:
        #prompt user for a score cutoff:
        input3 = int(input('enter the score cutoff you would like:\n '))
        function3(input3,cleanlist)
        
if __name__ == "__main__":
    main()
