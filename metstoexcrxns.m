function uniquemetstorxnindsornames = metstoexcrxns(metsarray,model,inds)
uniquemetstorxnindsornames=containers.Map;
for i=1:length(metsarray)
    if(strcmp(metsarray{i},'34hpp'))
        excrxnname='EX_34hpp';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    elseif(strcmp(metsarray{i},'glc_D'))
        excrxnname='EX_glc(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    elseif(strcmp(metsarray{i},'udpgal/udpg'))
        excrxnname1='EX_udpgal(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_udpg(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind1(1) excrxnind2(1)];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname1 excrxnname2};
        end
    elseif(strcmp(metsarray{i},'lac_D/lac_L'))
        excrxnname1='EX_lac_D(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_lac_L(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind1(1) excrxnind2(1)];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname1 excrxnname2};
        end
    elseif(strcmp(metsarray{i},'sbt_D'))
        excrxnname='EX_sbt-d(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    elseif(strcmp(metsarray{i},'tyr_l'))
        excrxnname='EX_tyr_L(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    elseif(strcmp(metsarray{i},'cit/icit'))
        excrxnname='EX_cit(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    elseif(strcmp(metsarray{i},'N/A'))
        excrxnname='';
        excrxnind=0;
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    else
        excrxnname=strcat('EX_',strcat(metsarray{i},'(e)'));
        excrxnind=find(strcmp(excrxnname,model.rxns));
        if(inds==1)
            uniquemetstorxnindsornames(metsarray{i})=[excrxnind];
        else
            uniquemetstorxnindsornames(metsarray{i})={excrxnname};
        end
    end
end
end

