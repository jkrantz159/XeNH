function parsave(alpha,beta,eta,Xed,Nd,Resfrac,LVfrac)

fid=fopen('success.txt','a+');
fprintf(fid,'%E,%E,%E,%E,%E,%E,%E',alpha,beta,eta,Xed,Nd,Resfrac,LVfrac);
fprintf(fid,'\n');
fclose(fid);