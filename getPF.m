for n=7
        for m=15
            parentPath='MaF';
            path=[int2str(m),'d'];
            mkdir(parentPath,path);
            problem=str2func(['MaF',int2str(n)]);
            Global=GLOBAL('-algorithm',{@KnEA,0.4},'-problem',problem,'-M',m);
            pf=Global.problem.PF(10000);
            name=[parentPath,'/',path,'/','MaF',int2str(n),'.pf'];
            writematrix(pf,name,'FileType','text','Delimiter','space')
        end
end