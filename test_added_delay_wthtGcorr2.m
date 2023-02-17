
%%%  //////////////////  Code for  Inflexion Point (down going)  \\\ MOH RAFIK /////////////////////////
%%%  ********************************************************************************************* %%%%
N=10;
mat =(0); outloop_count =0;
for p =2:N
    outloop_count = outloop_count+1 ;
    file = ['charge_',num2str(p),'.csv'];
    a = dlmread(file);
    [row,col] = size(a);
    a1=a(:,1); a2=a(:,2);

    a_min = min(a(:,2));
    
        h=((2.165 - 2.12)/((2E+4)-1));
    
    Ts = a1(2,1)-a1(1,1);
    Fs = (1/Ts);
    m = 30;    % window_size = ((m-1)*Ts)
    fprintf('window_size = %d\n',((m-1)*Ts));
    n = ceil(5e-9/((m-1)*Ts)); %  peak time is 20ns so we take 1/4 of peak time.((1/4)*(20ns))
    % n = ceil(10e-9/((m-1)*Ts)); %peak time is 20ns so we take 1/2 of peak time.((1/4)*(20ns))
    fprintf('n = %d\n',n);
    % m = 30; n = 10;
    % disp((m-1)*Ts);
    disp(n*(m-1)*Ts);
    infl_arr = zeros(m,1);
    min_infl_arr = zeros(n,1);
    %infxn_arr = zeros(N,1);
    
    % just before the break put the expresion for  inflexion point
    %   infx_point =  count*m - (m-1); where m is the size of window, means howmany points
    %   taken to be consideration. k means how many times  the loop performed.
    k=0;
    count=0;
    flag=0;
    for i= 1:m:row
        count = count+1;
        last_term = 1 +((count-1)*(m));
        if (last_term > row)
            break
        end
        q = i+m-1;
        if (q > row)
            q = row;
        end
        % infl_arr = a2(i:i+m-1,1);
        infl_arr = a2(i:q,1);
        k = k+1;
        
        %   disp(infl_arr);
        %     disp(min(infl_arr));
        %disp(k);
        if (k<=n)
            min_infl_arr(k,1)= min(infl_arr);
        end
        if k == n
            g=1;
            for j=2:n
                if ( min_infl_arr(j-1,1)> min_infl_arr(j,1))
                    g = g+1;
                    if ((g == n)&&(j == n))
                        %                     disp(min_infl_arr(j-2,1));
                        %                     disp(min_infl_arr(j-1,1));
                        disp(g);
                        disp(j);
                        %                     disp(flag);
                        flag = 1;
                        break
                    end
                end
            end
            % fprintf('just after  the if cond k = %d\t',k);
            k=0;
            %fprintf('just after the if cond after set k=0 k = %d\t\n',k);
        end
        if (flag == 1)
            infx_point = ((count.*m)-(m-1));
            %disp(infx_point);
            break
        end
    end
    mat = [mat,infx_point]; %#ok<AGROW>
    %disp(mat);
    % mat = mat(1,2:N+1);
    %x=a(mat(1,p+1),1);
    %x = a(infx_point,1);
    
    length = (a_max -a_min)/(h);
    length = ceil(length);
    x1_pconst = zeros(length,1);
    for i=1:length
        x1_pconst(i,1)= a(infx_point-(m*n),1);
        % x1_pconst(i,1)= a(infx_point,1);
    end
    % y1_pt1 = a_min: h : a_max;
    % figure
    % plot(a1,a2)
    % hold on
    % plot(x1_pconst,y1_pt1);
    % grid on
end
mat = mat(1,2:end);
mat = [mat(1),mat];
Down_inflx_point_unfltr = mat-(m*n);
Down_inflx_point_min_unfltrd = min(Down_inflx_point_unfltr);

%%%  ////////////////// End of Down inflexion point //////////////
%%%  **************************************************************************** %%%%


%%%  ///////////// applying and saving all the results hist,after and before standard deviation //////////////
close all
clc
fid = fopen('final1.txt','wt');
fprintf(fid,'\n');
fclose(fid);

N=10;
for p=1:N\
    filename = ['charge_',num2str(p),'.csv'];
    a = dlmread(filename);
    [row,col]=size(a);
    %   //* a is the data file which has two columns first column represents the time(ns) in 50ns scale // MOH RAFIK //
    %    and 2nd column represents the amplitude(Volt) of the output from from hardroc. // MOH RAFIK //   *//
    step_size = ((a(row,1)- a(1,1))/(row-1));  % time scale step size from the given data vector.
    % a1= a(:,1);   a2 = a(:,2);
    
    %//////// *********** index_min value of unfiltered signal ****************\\\\
    a1= a(:,1);  a2 = a(:,2);
    a_max = max(a(:,2));
    a_min = min(a(:,2));
    [val,idx] = min(a2);
    indx_minvalue_unfltrd = idx;
    peak_time_t_unfltrd = a(indx_minvalue_unfltrd,1);  % peak time of unfiltered signal.
    %        num_indx_dealy = ((peak_time_t_unfltrd/step_size)+1);
    % ***********  ////////////// min value index calculation is done \\\\\\\\\\\\ ************.
    
    
    % YMIN = min(a2); YMAX= max(a2);   %#ok<NASGU>
    for i=1:row
        if (a1(i)>= -1.5e-7)
            indx1 = i;
            break
        end
    end
    %  disp(indx1);
    %     i = 0;
    for i = row:-1:1
        if (a1(i) <= 1.5e-7)
            indx2 = i;
            break
        end
    end
    % disp(indx2);
    %a1 = a1(indx1:indx2);
    
    XMIN = a1(indx1); XMAX = a1(indx2);
    YMIN = min(a2(indx1:row)); YMAX= max(a2(indx1:row));
    
    [B,A] = butter(5,0.002,'low');  % here 5 is the order and 0.002 is the cuttoff freq.
    [h,w] = freqz(B,A,2001);
    
    y_butter_filtr = filter(B,A,a2);      % filtered data.
    y1_butter_filtr = y_butter_filtr(indx1:indx2);  % Y1 filtered data with /// TIME Axis limit /// the range.
    a1_y1_butter_filtr = a1(indx1:indx2);
    y1_butter_filtr_MAX = max(y1_butter_filtr(:,1));  % maximum of filtered data.
    [row_filt,col_filt]= size(y1_butter_filtr);
    y1_butter_filtr_MIN = min(y1_butter_filtr(:,1));  % minimum of filtered data.
    a_y_butter_filtr = [a1_y1_butter_filtr,y1_butter_filtr] ;
    [row_filt_vect,col_filt_vect] = size(a_y_butter_filtr);
    [val,idx] = min(y1_butter_filtr(:,1));
    indx_minvalue = idx;
    peak_time_t = a1_y1_butter_filtr(indx_minvalue,1);
    delay_filter_time = (peak_time_t - peak_time_t_unfltrd); % delay b/w real(Raw) and filtered pulse.
    indx_num_dealy = ((delay_filter_time/step_size)+1);
    disp('delay_filter_time =');
    disp(delay_filter_time);
    
    for i=1:row
        if ((a_y_butter_filtr(i,1)>= -0.025E-7)&&(a_y_butter_filtr(i,1)<= 0))
            halfend_offset = i;
            disp('i = ');
            disp(halfend_offset);
            break
        end
    end
    vect_mean = a_y_butter_filtr((1:halfend_offset),2);
    mean_offset = mean(vect_mean);
 
%end   
    Fs = 40e9;          % sampling frequency
    noise_after_filtering = std(y_butter_filtr(indx1+indx_num_dealy-1:Down_inflx_point_min_unfltrd+indx_num_dealy-1,1)); % standard deviation.
    %fprintf('sigma_after_filter = %e\n',noise_after_filtering);
    y_filter_peak2peak = max(y_butter_filtr(indx1+indx_num_dealy-1:Down_inflx_point_min_unfltrd+indx_num_dealy-1))-min(y_butter_filtr(indx1+indx_num_dealy-1:Down_inflx_point_min_unfltrd+indx_num_dealy-1));
    noise_after_filtering = noise_after_filtering*1E+3;
    noise_original_signal = std(a2(indx1:Down_inflx_point_min_unfltrd,1)); % standard deviation.
    
	
	
    y_rawsignal_peak2peak = max(a2(indx1:Down_inflx_point_min_unfltrd,1))- min(a2(indx1:Down_inflx_point_min_unfltrd,1));
    %fprintf('sigma_orignal_signal = %e\n',noise_original_signal);
    noise_original_signal = noise_original_signal*1E+3 ;
    noise_improved = (noise_original_signal - noise_after_filtering);
    %fprintf('sigma_after_filter = %e\n',noise_improved);
    improvement_signal_peak2peak =  (y_rawsignal_peak2peak - y_filter_peak2peak);
    
    figure
    ax1 = subplot(1,2,1);
    %hist(y1_butter_filtr(indx1:Down_inflx_point_min),500);
    hist(y1_butter_filtr(indx1+indx_num_dealy-1:Down_inflx_point_min_unfltrd+indx_num_dealy-1),100);
    annotation('textbox',[0.10 0.25 0.1 0.10],'String',{['Sigma =' num2str(noise_after_filtering)]});
    %% annotation('textbox',[x y w h]) creates an editable text box annotation with its lower left corner at the point x,y,
    %% a width w, and a height h, specified in normalized figure units. Specify x, y, w, and h in a single vector.
    title('After filtering')
    ax2 = subplot(1,2,2);
    %hist(a2(indx1:Down_inflx_point_min),500);
    hist(a2(indx1:Down_inflx_point_min_unfltrd),100);
    annotation('textbox',[0.60 0.65 0.1 0.10],'String',{['Sigma =' num2str(noise_original_signal)]});
    title('Before filtering')
    linkaxes([ax2,ax1],'xy');
    set(gcf,'Units','inches');
    set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
    screenposition = get(gcf,'Position');
    %%%set(gca,'FontSize',10.2,'fontWeight','bold');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    % title('VOLT(KV) Vs EFFICIENCY(%)'); %title('DAC2 Vs Qinj( pC )');
    %  xlabel( Voltage(KV)');
    %  xlim([9.2 11.6]);
    %print -dpdf -painters sigma_before_after
    % print(gcf,'-depsc','sigma_before1_after.eps');
    print( '-dpdf','-painters',strcat('sigma_before_after',num2str(p)));
    % print(gcf,'-depsc','only_crosstalk_mainstripsf61.eps');
    
    figure
    % subplot(1,2)
    ax1 = subplot(2,1,1);
    %hist(y1_butter_filtr(indx1:Down_inflx_point_min),500);
    hist(y1_butter_filtr(indx1+indx_num_dealy-1:Down_inflx_point_min_unfltrd+indx_num_dealy-1),100);
    title('After filtering')
    ax2 = subplot(2,1,2);
    %hist(a2(indx1:Down_inflx_point_min),500);
    hist(a2(indx1:Down_inflx_point_min_unfltrd),100);
    title(' Before filtering')
    linkaxes([ax2,ax1],'xy');
    set(gcf,'Units','inches');
    set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
    screenposition = get(gcf,'Position');
    %%set(gca,'FontSize',10.2,'fontWeight','bold');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    % title('VOLT(KV) Vs EFFICIENCY(%)'); %title('DAC2 Vs Qinj( pC )');
    %  xlabel( Voltage(KV)');
    %    xlim([9.2 11.6]);
    %print -dpdf -painters sigma_before_aftercolumn
    % print(gcf,'-depsc','sigma_before2_after_coulmn.eps');
    print( '-dpdf','-painters',strcat('sigma_before_aftercolumn',num2str(p)));
    
    figure
    subplot(2,1,1);
    plot(a1(indx1:indx2),y1_butter_filtr,'k');
    %     xlim([-1.5E-7 1.5E-7])
    %     ylim([min(y1_butter_filtr) max(y1_butter_filtr)])
    grid on
    subplot(2,1,2);
    plot(a1(indx1:indx2),a2(indx1:indx2));
    %     xlim([-1.5E-7 1.5E-7])
    %     ylim([min(a2(indx1:indx2)) max(a2(indx1:indx2))])
    grid on
    set(gcf,'Units','inches');
    % set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
    screenposition = get(gcf,'Position');
    %%set(gca,'FontSize',10.2,'fontWeight','bold');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    % title('VOLT(KV) Vs EFFICIENCY(%)'); %title('DAC2 Vs Qinj( pC )');
    %  xlabel( Voltage(KV)');
    %    xlim([9.2 11.6]);
    %  print -dpdf -painters noise_noise_free_signal
    %    print(gcf,'-depsc','noise_free_signal.eps');
    print( '-dpdf','-painters',strcat('noise_noise_free_signal',num2str(p)));
    
    figure
    grid on
    plot(a1(indx1:indx2),y1_butter_filtr,'ok','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',1.5);
    grid on
    %     xlim([-1.5E-7 1.5E-7])
    %     ylim([min(y1_butter_filtr) max(y1_butter_filtr)])
    hold on
    plot(a1(indx1:indx2),a2(indx1:indx2),'*g','LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor',[0 0 0],'MarkerSize',1.5);
    set(gcf,'Units','inches');
    % set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
    screenposition = get(gcf,'Position');
    %%set(gca,'FontSize',10.2,'fontWeight','bold');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    % title('VOLT(KV) Vs EFFICIENCY(%)'); %title('DAC2 Vs Qinj( pC )');
    %  xlabel( Voltage(KV)');
    %    xlim([9.2 11.6]);
    % print -dpdf -painters superimpose_noise_free_signal
    %print(gcf,'-depsc','superimpose_noise_free_signal.eps');
    grid on
    print(gcf,'-depsc',strcat('superimpose_noise_free_signal',num2str(p)));
    print( '-dpdf','-painters',strcat('superimpose_noise_free_signal',num2str(p)));
    
    %end         here end is commented.
    
    %     fprintf('halfend_offset = %d\n',halfend_offset);
    %     fprintf('mean_offset = %d\n',mean_offset);
    
    l = (mean_offset - y1_butter_filtr_MIN);
    y2_butter_filtr = (mean_offset - 0.9*l);  % Amplitude at 90% of the max value. it was y2 earlier.  y2_butter_filtr
    y3_butter_filtr = (mean_offset - 0.1*l);  % Amplitude at 10% of the max value. it was y1 earlier.
    y4_butter_filtr = (mean_offset - 0.05*l); %  For peaking time. Amplitude at 5% of the max value.
    
    % /////// For the first_FWHM calculation ///// MOH RAFIK /////
    
    y_FWHM = ((mean_offset + y1_butter_filtr_MIN)/2);
    y1_butter_filtr_FWHM =(y1_butter_filtr - y_FWHM);
    y1_butter_filtr_FWHMpos = abs(y1_butter_filtr - y_FWHM);
    y1_butter_filtr_FWHMpos_untilminvalindx = y1_butter_filtr_FWHMpos(1:indx_minvalue,1);
    [val,idx]= min(y1_butter_filtr_FWHMpos_untilminvalindx);
    indx_t_width1 = idx;   %this will get corresponding time value.% maximum of filtered data.
    [row_filt,col_filt]= size(y1_butter_filtr);
    y1_butter_filtr_MIN = min(y1_butter_filtr);
    min1_y1_butter_filtr_FWHMpos = val;
    t_width1 = a1_y1_butter_filtr(indx_t_width1,1);
    
    
    % /////// For the second_FWHM after the indx_minvalue /////
    
    y1_butter_filtr_FWHMpos_min2lastelmnt_vect = y1_butter_filtr_FWHMpos((indx_minvalue:row_filt),1);
    [val,idx] = min(y1_butter_filtr_FWHMpos_min2lastelmnt_vect);
    indx_t_width2 = indx_minvalue + idx-1;
    min2_y1_butter_filtr_FWHMpos = val;
    t_width2 = a1_y1_butter_filtr(indx_t_width2,1);
    pulse_width = (t_width2 - t_width1);
    %  ////// end of FWHM calculation ///////
    
    %***************************************************************************************///
    % ///// NEW rise time 10p calculation.////
    Vect_rise_or_fall_10p = ( y1_butter_filtr - y3_butter_filtr );
    for i= indx_minvalue:-1:1
        if(i<indx_t_width1)
            if (Vect_rise_or_fall_10p(i,1)>=0)
                v = i;
                break
            end
        end
    end
    %fprintf('v = %d\n',v);
    rise10p_vect = zeros(7,2);
    for i= v-6:v
        for j=1:2
            rise10p_vect(i-v+7,j) = a_y_butter_filtr(i,j);
        end
    end
    [val,idx] = min(abs(rise10p_vect(1:7,2)));
    idx_rise_time10p = idx-1 + v-6 ;
    t1_d  =  a_y_butter_filtr(idx_rise_time10p,1);
    %fprintf('t1_d = %e\n',t1_d);
    % % ///// End of NEW rise time 10p calculation.////
    
    %***************************************************************************************///
    % ///// New Fall time 10p calculation ////
    Vect_rise_or_fall_10p = (y1_butter_filtr - y3_butter_filtr);
    for i= indx_minvalue:indx2
        if(i>indx_t_width2)
            if (Vect_rise_or_fall_10p(i,1)>=0)
                v_fall_10p = i;
                break
            end
        end
    end
    % fprintf('v_fall_10p = %d\n',v_fall_10p);
    fall10p_vect = zeros(7,2);
    for i = v_fall_10p:v_fall_10p+6
        for j = 1:2
            fall10p_vect(i-v_fall_10p+1,j) = a_y_butter_filtr(i,j);
        end
    end
    [val,idx] = min(abs(fall10p_vect(1:7,2)));               %#ok<*ASGLU>
    idx_fall_time10p = idx + v_fall_10p-1;
    Fall_time10p_t2  =  a_y_butter_filtr(idx_fall_time10p,1);
    %fprintf('Fall_time10p_t2 = %e\n',t1_d);
    % /////End of NEW Fall time 10p calculation.////
    %***************************************************************************************///
    
    
    %****************************************************************************************************************%
    % ////  This following block for the RISE time 90p calculation. ////
    % %     switch p
    % %         case 1
    Vect_rise_or_fall_90p = (y1_butter_filtr - y2_butter_filtr);
    vector_for90p_min = abs(y1_butter_filtr- y2_butter_filtr);
    Vect_rise_90p_selctd = vector_for90p_min((indx_t_width1:indx_minvalue),1);
    [val,idx] = min(Vect_rise_90p_selctd);
    y2_zero = val;
    indx_rise_90p_risingedge = idx;
    indx_rise_90p = indx_t_width1 + indx_rise_90p_risingedge-1; %% indx_rise_90p will represent the index corresponding to the 90% for rise time.
    t2_d =a_y_butter_filtr(indx_rise_90p,1);  % Time is directly taken from first col of the data according to the 90% value.
    
    %     if p == 3
    %         fprintf('indx_rise_90p = %e\n',indx_rise_90p);
    %         fprintf('t2_d = %e\n',t2_d);
    %         fprintf('t1_d = %e\n',t1_d);
    %     end
    %\\\\\\\\\ End OF block RISE TIME 90P //////////
    %******************************************************************************************************%
    
    % /////  Fall time 90p calculation BLOCK ////
    vect_Fall90p = vector_for90p_min((indx_minvalue:indx_t_width2),1);
    [val,idx] = min(vect_Fall90p);
    indx_fall_90p = indx_minvalue + idx -1 ;
    Fall_time90p_t1 = a_y_butter_filtr(indx_fall_90p,1);
    Fall_time = (Fall_time10p_t2 - Fall_time90p_t1) ;
    
    %\\\\\\\\\ End OF block Fall TIME 90P //////////
    %******************************************************************************************************%
    
    %     t2 = -250E-9 + (500E-9/(indx2-1))*(indx_rise_90p-1);
    %     t1 = -250E-9 + (500E-9/(indx2-1))*(idx1_of_corres_a - 1);
    %     Rise_time = (t2-t1);
    Rise_time_d = (t2_d - t1_d);   % Rise time of a given pulse.According to the given time directly from the data files without using own formula.
    % fprintf('Rise_time_d = %e\n',Rise_time_d);
    
    %******************************************************************************************************%
    
    % ///// Peaking_time ////   %Peaking_time = (peak_time_t) - (t1_d);
    
    peaking_time_vect = (y1_butter_filtr - y4_butter_filtr );   % indx_minvalue
    for i= indx_minvalue:-1:1
        if(i<indx_t_width1)
            if (peaking_time_vect(i,1)>=0)
                v_peak_value = i;
                break
            end
        end
    end
    peak_time_vect = zeros(7,2);
    for i= v-6:v
        for j=1:2
            peak_time_vect(i-v+7,j) = a_y_butter_filtr(i,j);
        end
    end
    [val,idx] = min(abs(peak_time_vect(1:7,2)));
    idx_peak_time_vect = idx-1 + v-6 ;
    peak_time_5p = a_y_butter_filtr(idx_peak_time_vect,1);
    % size_peaking_time_vect = 200;
    % peaking_time_vect = zeros(size_peaking_time_vect,2);
    %     for i = idx_rise_time10p-(size_peaking_time_vect-1):1:idx_rise_time10p
    %         for j = 1:2
    %             peaking_time_vect(i- idx_rise_time10p + size_peaking_time_vect,j) = a_y_butter_filtr(i,j);
    %         end
    %     end
    %     [val,idx] = min(abs(peaking_time_vect(1:size_peaking_time_vect,2)-(mean_offset)));
    %    peak_firstpoint = (idx -1) + (idx_rise_time10p - size_peaking_time_vect-1);
    %         peak_time_firstpoint = a_y_butter_filtr(peak_firstpoint,1);
    Peaking_time = (peak_time_t) - (peak_time_5p);
    
    %\\\\\\\\\ End OF Peaking time BLOCK //////////
    %******************************************************************************************************%
    
    h=((2.165 - 2.12)/((2E+4)-1));
    
    % % For appending the data continously after extracting the data from DSO each files.
    fid = fopen('final1.txt','at');
    fmt = '%e %e %e %e %e %e\n';
    A = [Rise_time_d Peaking_time Fall_time pulse_width delay_filter_time indx_num_dealy];
    fprintf(fid,fmt,A);
    fclose(fid);
    
    y1_p=zeros(indx2,1);
    for i=1:indx2
        y1_p(i,1)= y3_butter_filtr;             % y3_butter_filtr = 10% of  its  max value. horizontal  line for 10%.
    end
    y2_p=zeros(indx2,1);
    for i=1:indx2
        y2_p(i,1)= y2_butter_filtr;               %y2 = 90% of  its  max value. horizontal  line for 90%.
    end
    length = (y1_butter_filtr_MAX -y1_butter_filtr_MIN)/(h);
    length = ceil(length);
    x1_pconst = zeros(length,1);
    for i=1:length
        x1_pconst(i,1)= t1_d;
    end
    y1_pt1 = y1_butter_filtr_MIN: h : y1_butter_filtr_MAX;   % y1_pt1 represnts vertical line at 10% of its max value at rising edge .
    
    x2_pconst = zeros(length,1);
    for i=1:length
        x2_pconst(i,1)= t2_d;
    end
    y2_pt1 = y1_butter_filtr_MIN: h : y1_butter_filtr_MAX;     % y2_pt1 represnts vertical line at 90% of its max value at rising edge .
    
    x3_pconst = zeros(length,1);
    for i=1:length
        x3_pconst(i,1)= t_width1;
    end
    y3_pt1width = y1_butter_filtr_MIN: h : y1_butter_filtr_MAX;   % y3_pt1width vertical line at the point where the FWHM exits at the rising edge.
    
    x4_pconst = zeros(length,1);
    for i=1:length
        x4_pconst(i,1)= t_width2;
    end
    y4_pt2width = y1_butter_filtr_MIN: h :y1_butter_filtr_MAX;    % y4_pt1width vertical line at the point where the FWHM exits at the Falling edge.
    
    
    x5_pconst = zeros(length,1);
    for i=1:length
        x5_pconst(i,1)= Fall_time90p_t1;
    end
    y5_FallTime_90p = y1_butter_filtr_MIN: h :y1_butter_filtr_MAX;     % y5_FallTime_90p represnts vertical line at 10% of its max value at rising edge .
    
    
    x6_pconst = zeros(length,1);
    for i=1:length
        x6_pconst(i,1)= Fall_time10p_t2;
    end
    y6_FallTime_10p = y1_butter_filtr_MIN: h :y1_butter_filtr_MAX;         % y5_FallTime_10p represnts vertical line at 10% of its max value at rising edge .
    
    
    x7_pconst = zeros(length,1);
    for i=1:length
        x7_pconst(i,1)= peak_time_5p;
    end
    y7_FallTime_10p = y1_butter_filtr_MIN: h :y1_butter_filtr_MAX;
    
    %*******************************************************************************************************************************////
    % ////// For the Original data plot (Pulse of OUT_Fsb).\\\\\\\
    figure
    % x1 = (1E-9)*(a1);
    plot(a1_y1_butter_filtr,y1_butter_filtr,'-ob','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor',[0 0 0],'MarkerSize',1.2);
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    % set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
    set(gca,'FontSize',10.2,'fontWeight','bold')
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    title('Time(ns) Vs Amplitude(mV)');
    xlabel('Time(ns)');
    xlim([-150E-9 150E-9]);
    axis tight
    ylabel('Amplitude(V)');
    % ylim([y1_butter_filtr_MIN y1_butter_filtr_MAX]);
    %   ylim([1.9 2.25]);
    %*******************************************************************************************************************************////
    % //////  plot of horizontal line of 10%.\\\\\\\
    hold on;
    plot(a1(indx1:indx2),y1_p(indx1:indx2),'r');                 %horizontal  line for 10%.
    
    %*******************************************************************************************************************************////
    % //////  plot of horizontal line of 90%.\\\\\\\
    hold on;
    plot(a1(indx1:indx2),y2_p(indx1:indx2),'r');               %horizontal  line for 90%.
    
    %*******************************************************************************************************************************////
    %//// y1_pt1 represnts vertical line at 10% of its max value at rising edge////.
    hold on;
    plot(x1_pconst, y1_pt1,'r', 'LineWidth', 0.5) ;
    %*******************************************************************************************************************************////
    %//// y2_pt1 represnts vertical line at 90% of its max value at rising edge////.
    hold on;
    plot(x2_pconst, y2_pt1, 'r','LineWidth', 0.5) ;
    %*******************************************************************************************************************************////
    %//// % Vertical line on rising edge for the pulse width ////.
    hold on;
    plot(x3_pconst, y3_pt1width,'k', 'LineWidth', 0.5) ;    % Vertical line on rising edge for the pulse width.
    
    %*******************************************************************************************************************************////
    %//// % Vertical line on falling edge for the pulse width ////.
    
    hold on;
    plot(x4_pconst, y4_pt2width, 'k','LineWidth', 0.5) ;    % Vertical line on falling edge for the pulse width.
    
    %*******************************************************************************************************************************////
    %//// % Vertical line at FallTime_90p ////.
    hold on;
    plot(x5_pconst, y5_FallTime_90p,'g', 'LineWidth', 0.5) ;
    
    %*******************************************************************************************************************************////
    %////  Vertical line at FallTime_10p ////.
    hold on;
    plot(x6_pconst, y6_FallTime_10p, 'g','LineWidth', 0.5) ;
    %hold on;
    grid on;
    
    
    hold on;
    plot(x7_pconst, y7_FallTime_10p, 'k','LineWidth', 0.5) ;
    %hold on;
    grid on;
    
    
    % hLine = imline(gca,[-2.5E-7 2.12],[ 2.5E-7 2.165]);
    
    %    print( '-dpdf','-painters',strcat('charge',num2str(p)));
    % print(gcf,'-depsc','only_crosstalk_mainstripsf61.eps');
    
end
q = dlmread('final1.txt');
y = zeros(10,6);
x=[ 10
    20
    30
    40
    50
    60
    70
    80
    90
    100
    ];

for i=1:10
    for j=1:6
        y(i,j)= (q(i,j))*1E+9;
    end
end
%disp(y);
y1 = y(:,1);   % Rise_time_d
y2 = y(:,2);   % Peaking_time
y3 = y(:,3);   % Fall_time
y4 = y(:,4);   % pulse_width
y5 = y(:,5);   % Filtered_delay
mean_delay_filtr = mean(y5);
figure
plot(x,y1,'-^r','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',4);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gca,'XMinorTick','on','YMinorTick','on');        % for minor ticks
set(gca,'FontSize',10.2,'fontWeight','bold')
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
% title('VOLT(KV) Vs EFFICIENCY(%)');
xlabel('Charge(fC)');
xlim([10 100]);
ylabel('Time(ns)');
ylim([5 35]);
grid on
%grid minor
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% ///////////////////// -------1
hold on;
plot(x,y2,'-ok','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',4);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
set(gca,'FontSize',10.2,'fontWeight','bold')
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
% title('VOLT(KV) Vs EFFICIENCY(%)');
% xlabel('EFFECTIVE VOLTAGE (KV)');
%  xlim([9.6 11.6]);
grid on
hold on;
plot(x,y3,'-<b','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',4);
set(gcf,'Units','inches');
set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
xlabel('Charge(fC)');
grid on
hold on;

plot(x,y4,'-dg','LineWidth',2.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',4);
set(gcf,'Units','inches');
set(gca,'XMinorTick','on','YMinorTick','on');  % for minor ticks
screenposition = get(gcf,'Position');
set(gca,'FontSize',10.2,'fontWeight','bold');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
% title('VOLT(KV) Vs EFFICIENCY(%)'); %title('DAC2 Vs Qinj( pC )');
xlabel('Charge(fC)');
grid on
hold on
plot(x,y5,'-*','LineWidth',2.5,'MarkerEdgeColor','c','MarkerFaceColor',[0 0 0],'MarkerSize',4);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gca,'XMinorTick','on','YMinorTick','on');        % for minor ticks
set(gca,'FontSize',10.2,'fontWeight','bold')
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
% title('VOLT(KV) Vs EFFICIENCY(%)');
annotation('textbox',[0.85 0.85 0.1 0.10],'String',{['Mean filter delay =' num2str(mean_delay_filtr)]});
xlabel('Charge(fC)');
xlim([10 100]);
ylabel('Time(ns)');
ylim([5 35]);
grid on
legend('Rise time','Peaking time','Fall time','Pulse width','Filtered delay','Location','best');
print -dpdf -painters RisepeakFallWidth_withNoise_Fsb0
print(gcf,'-depsc','RisepeakFallWidth_withNoise_fsb0.eps')






