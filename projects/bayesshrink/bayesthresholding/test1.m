%test.m
clc;clear;close all;

m=30;
fb=200;
p=5;
q=20;

t=linspace(-0.2,0.2,200);
y=1/(q-p)*sqrt(fb) *(sinc(fb*t/m)).^m.*(q*sinc(2*q*t)-p*sinc(2*p*t));

figure;plot(t,y);