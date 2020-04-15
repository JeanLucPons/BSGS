# Baby Step Giant Step for SECPK1

A simple Baby Step Giant Step program for SecpK1.

Structure of the input file:
* All values are in hex format
* Public keys can be given either in compressed or uncompressed format

```
Babystep size
Start range
End range
Key #1
Key #2
...
```

ex

```
40000000
49dccfd96dc5df56487436f5a1b18c4f5d34f65ddb48cb5e0000000000000000
49dccfd96dc5df56487436f5a1b18c4f5d34f65ddb48cb5effffffffffffffff
0459A3BFDAD718C9D3FAC7C187F1139F0815AC5D923910D516E186AFDA28B221DC994327554CED887AAE5D211A2407CDD025CFC3779ECB9C9D7F2F1A1DDF3E9FF8
0335BB25364370D4DD14A9FC2B406D398C4B53C85BE58FCC7297BD34004602EBEC
```

# How it works

It uses a hash table to store the baby steps (the generator table).
2^30 baby steps is about 9GB (we do not store full point, only a part of it, false collisions are resolved later).
Here is a brief description of the algoritm:

We have to solve P = k.G, we know that k lies in the range ]k1,k2]

```
 m = Baby Step Size # With enough memory, the best choice is m = sqrt(k2-k1)
 Compute a table of [G,2G,..,b.G,...,m.G]
 S = P - k1.G
 found = false
 step = 0
 while not found and step<(k1-k2) {
   if S is in the table {
     k = k1 + step + b
     found = true
   }
   S = S - m.G
   step = step + m
 }
```








