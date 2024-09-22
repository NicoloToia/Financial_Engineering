# Q1

200 obligors, no correlation

``` ext
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (200 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2160; down prob. (input) 0.2160 
stay prob. (simulated) 0.7791; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.000 (200 names)–––
VaR - default only 1.49 
VaR - default and migration 1.95 
```

# Q2

## Case 1

200 obligors, correlation $\sqrt{0.12}$

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (200 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2162; down prob. (input) 0.2160 
stay prob. (simulated) 0.7788; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.120 (200 names)–––
VaR - default only 3.58 
VaR - default and migration 4.68 
```

## Case 2

200 obligors, correlation $\sqrt{0.24}$

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (200 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2164; down prob. (input) 0.2160 
stay prob. (simulated) 0.7786; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.240 (200 names)–––
VaR - default only 6.85 
VaR - default and migration 8.15 
```

## Case 3

200 obligors, correlation with IRB formula

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (200 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2163; down prob. (input) 0.2160 
stay prob. (simulated) 0.7787; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.213 (200 names)–––
VaR - default only 6.26 
VaR - default and migration 7.48 
```

# Q3

## Case 1

Use 20 obligors, correlation $\sqrt{0.12}$

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (20 names)–––
PD         (simulated) 0.0049; PD         (input) 0.0050
down prob. (simulated) 0.2160; down prob. (input) 0.2160 
stay prob. (simulated) 0.7791; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.120 (20 names)–––
VaR - default only 5.96 
VaR - default and migration 7.31 
```

## Case 2

Use 20 obligors, correlation $\sqrt{0.24}$

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (20 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2164; down prob. (input) 0.2160 
stay prob. (simulated) 0.7786; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.240 (20 names)–––
VaR - default only 8.94 
VaR - default and migration 10.39 
```

## Case 3

Use 20 obligors, correlation with IRB formula

```text
––– Part I Q1: Present Value in a years’time of the IG bond –––
1y fwds: 100.51 IG; 97.66 HY; 40.00 Def 
Present Value in a years time: 99.59
 
––– Test of the accuracy of the MC simulation (20 names)–––
PD         (simulated) 0.0050; PD         (input) 0.0050
down prob. (simulated) 0.2163; down prob. (input) 0.2160 
stay prob. (simulated) 0.7787; stay prob. (input) 0.7790 
up prob.   (simulated) 0.0000; up prob.   (input) 0.0000 
 
––– Part I Q2/3: Credit VaR with AVR correlation 0.213 (20 names)–––
VaR - default only 8.94 
VaR - default and migration 10.20 
```