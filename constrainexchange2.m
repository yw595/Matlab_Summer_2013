function model = constrainexchange2(model1,cellline)
inputFI=fopen('../Recon2 Constraints.csv','r');
excnumarrayconstraints=[];

line=fgetl(inputFI);
i=0;
while line~=-1
    words=strsplit(line,',');
    i=i+1;
    if(i>=11 && i<=101)
        excnumarrayconstraintsrow=[];
        for j=9:2:127
            %if(j~=11 && j~=89)
                excnumarrayconstraintsrow=[excnumarrayconstraintsrow str2num(words{j})];
            %end
        end
        excnumarrayconstraints(end+1,:)=excnumarrayconstraintsrow;
    end
    line=fgetl(inputFI);
end

[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
uniquemetstoexcrxnnames=metstoexcrxns(metsarray,model1,2);
celllinesarray=exctextarray(9,10:2:128);

exctokeep={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
celllineind=find(strcmp(cellline,celllinesarray));
model=model1;
excsystemindices=find(ismember(model.subSystems,'Exchange/demand reaction'));
for i=1:length(excsystemindices)
    %model.lb(excsystemindices(i))=0;
    %model.ub(excsystemindices(i))=0;
    %disp(model.rxns{excsystemindices(i)});
end
for i=1:length(model.rxns)
    if (~isempty(regexp(model1.rxns{i},'^EX(.)*\(e\)$')))
        %if (sum(strcmp(model1.rxns{i},exctokeep))==0)
            model.lb(i)=0;
            %model.ub(i)=0;
        %end
    end
end
excnumarrayconstraintscol=excnumarrayconstraints(:,celllineind);
%excnumarrayconstraintscol
for j=1:length(excnumarrayconstraintscol)
    met=metsarray{j};
    %met
    excrxnname=uniquemetstoexcrxnnames(met);
    excrxnind=find(ismember(model.rxns,excrxnname));
    if(excnumarrayconstraintscol(j)>0)
        for k=1:length(excrxnind)
            %model.ub(excrxnind(k))=model.ub(excrxnind(k))+excnumarrayconstraintscol(j)/length(excrxnind);
        end
    elseif(excnumarrayconstraintscol(j)<0)
        for k=1:length(excrxnind)
            if(sum(find(ismember(exctokeep,model.rxns{excrxnind(k)})))~=0)
                model.lb(excrxnind(k))=model.lb(excrxnind(k))+excnumarrayconstraintscol(j)/length(excrxnind);
            end
        end
    end
end
for j=1:length(exctokeep)
    ind=find(ismember(model.rxns,exctokeep{j}));
    if(model.lb(ind)==0)
        model.lb(ind)=-1000;
    end
    disp([num2str(ind) ' ' num2str(model.rxns{ind}) ' ' num2str(model.lb(ind)) ' ' num2str(model.ub(ind))]);
end
for j=1:length(model.rxns)
    if (~isempty(regexp(model1.rxns{j},'^EX(.)*\(e\)$')))
        %disp([num2str(j) ' ' num2str(model.rxns{j}) ' ' num2str(model.lb(j)) ' ' num2str(model.ub(j))]);
    end
end
indstoshow=find(model.lb);%model.rxns(find(model.lb))
for i=1:length(indstoshow)
    if (~isempty(regexp(model1.rxns{indstoshow(i)},'^EX(.)*\(e\)$')))
        %disp([num2str(indstoshow(i)) ' ' num2str(model.rxns{indstoshow(i)}) ' ' num2str(model.lb(indstoshow(i))) ' ' num2str(model.ub(indstoshow(i)))]);
    end
end
%find(model.lb)
%model.lb(find(model.lb))
%model.rxns(find(model.ub))
%model.ub(find(model.ub))
end

