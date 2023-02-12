function [ out ] = Abbreviate_str2num( in )
cd('E:\Andrew system');
Converter_str2num(in); 
eval('out=simple_str2num(in);');
delete simple_str2num.m
end                                                                                                                                               