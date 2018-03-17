import sys, os, fcntl

flags = fcntl.fcntl(sys.stdout, fcntl.F_GETFL)
fcntl.fcntl(sys.stdout, fcntl.F_SETFL, flags&~os.O_NONBLOCK)

with open(sys.argv[1], 'r') as f:
    print f.read()
