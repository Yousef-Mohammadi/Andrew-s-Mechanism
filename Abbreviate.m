function [ out ] = Abbreviate( in )
cd('E:\Andrew system');
Converter(in); 
eval('out=simple(in);');
delete simple.m
end