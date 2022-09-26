Testing _Float16

Macos with M1/M2 CPU:
```
make -f Makefile.mac
sudo ./a.out
```

Linux with g++-12
```
make -f Makefile.linux
sudo ./a.out
```

ARMv8(AWS graviton3/Ubuntu 20.04LTS)
```
make -f Makefile.graviton3
sudo ./a.out
```
