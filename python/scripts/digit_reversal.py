import math

a = 123412840

# a = 123456789


b = 0
ndig = math.ceil(math.log10(a)) # Get the number of digits
for idx in range(ndig):
    b = 10*(b + 10*(a/10 - math.floor(a/10)))
    a = math.floor(a/10)

b = round(b/10)
print(b)