clear all;
clc;
Hname = 'H_336_672_nonReg.mat';
load(Hname)

shag = 0.1;
[SmartPattern, Rates] = GroupAndSort2(H,0.5,shag);
% SmartPattern = SmartPattern(2:end);
if (Rates(end)-Rates(end-1))<shag/2
   Rates = Rates(1:end-1)
end

for i = 1:length(SmartPattern)
   if Rates(i)==0.5
      SNR = 1:0.2:3
   elseif Rates(i)==0.6
      SNR = 2:0.2:3.8
   elseif Rates(i)==0.7
      SNR = 3:0.2:4.4
   elseif Rates(i)>0.7
      SNR = 3:0.2:6
   end
   
   BER = PLATFORM(H, Rates(i), SNR, SmartPattern{i});
   CurveName = strcat('BER_SMART_', num2str(Rates(i)),'_',Hname);
   SNRName = strcat('SNR_', num2str(Rates(i)),'_',Hname);
   save(CurveName,'BER');
   save(SNRName,'SNR');
   
   if i>1
      [SundromeBasedPattern, ~] = find_puncture_pattern(H, length(SmartPattern{i}), 1000);
      BER = PLATFORM(H, Rates(i), SNR, SundromeBasedPattern);
      CurveName = strcat('BER_IITP_', num2str(Rates(i)),'_',Hname);
      save(CurveName,'BER');
   end
end
      