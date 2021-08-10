x="abaasdasdusduhfikliilkjhgjhgjhgh"

def checkpalindrome(s,i):
    if len(s)>2:
        rev=s[::-1]
        if rev==s:
            print(i,":",s)
i=0
for l in x:
    s=""
    k=i
    while k < len(x):
        s=s+x[k]
        checkpalindrome(s,i)
        k=k+1
    i=i+1
