function plot_pantompkins_qrs_detection(testdata,filt_dat,int_dat,qrs_pos,thF1,thI1,fs)
    figure();
    L = 10000;
    t = (0:L-1) / fs;
    plot(t,testdata(1:L),'-k');
    hold on;
    qrs_pos = qrs_pos(qrs_pos<=L);
    if (~isempty(qrs_pos))
    plot((qrs_pos-1)/fs,testdata(qrs_pos),'og');
    end
    grid on;
    title("R peak detection result on Lead 1")
    %plot(filt_dat,'-m');
    %plot(int_dat,'-m');

    %plot(qrs_pos,thF1,'-m');
    %plot(qrs_pos,thI1,'-m');
    
    
    hold off;
end