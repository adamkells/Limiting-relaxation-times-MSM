%% Code by Adam Kells to fit to analytic data and compute errors
close all
clear all

for i=[1:10]
    load(['../HMM_final/new_data/data_3_' num2str(i) '.mat'])
        mm(i,:)=mm_all;
        hmm(i,:)=hmm_all;
        hmm3(i,:)=hmm3_all;
end

% remove the first data point since the fitting procedure is not valid in
% the very short lag time limit
mm_all=mm(:,2:end);
hmm_all=hmm(:,2:end);
hmm3_all=hmm3(:,2:end);
rel_exact=4802.3;

t=double(lags(2:end))*0.25;
for j=1:10
    [x1]=optimize_fit(t,mm_all(j,:)');
    p1(j)=x1(2);
end
mm_fit_err=std(p1); % error on the fit

% average relaxation time over all the runs
mm=mean(mm_all);
hmm2=mean(hmm_all);
hmm_3=mean(hmm3_all);

% standard deviation on the rel time over the runs
mm_err=std(mm_all);
hmm2_err=std(hmm_all);
hmm3_err=std(hmm3_all);

[x1]=optimize_fit(t,mm');

%% making figures with data calculated above
plot(t,rel_exact*ones(length(t),1),'color','k','linewidth',2)
hold on
% original data with errors
shadedErrorBar(t,mm,mm_err,'lineprops','b*')
shadedErrorBar(t,hmm2,hmm2_err,'lineprops','g*')

%keyboard
t_new=linspace(t(1),t(end),1000);
plot(t_new,(t_new.*mm_fit(2)./(t_new+mm_fit(2)*mm_fit(1))))

shadedErrorBar(t,mm_fit(2)*ones(length(t),1),max(mm_fit_err)*ones(length(t),1),'lineprops','r')%,'linewidth',2)
xlabel('\tau')
ylabel('\mu_2^{relax-MSM}')
h = findobj(gca,'Type','line');
legend([h(11) h(8) h(5) h(4) h(1)],{'Exact','MSM','HMM','Best fit','Limit of fit'})
shadedErrorBar(t,hmm_3,hmm3_err,'lineprops','k')%,'linewidth',2)

saveas(gcf,'analytic.fig','fig')
saveas(gcf,'analytic.png','png')