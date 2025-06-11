function plot_pantompkins_qrs_detection(testdata,filt_dat,int_dat,qrs_pos,thF1,thI1)
    figure();
    plot(testdata,'-b');
    hold on;
    if (~isempty(qrs_pos))
    plot(qrs_pos,testdata(qrs_pos),'og');
    end
    plot(filt_dat,'-m');
    %plot(int_dat,'-m');

    plot(qrs_pos,thF1,'-m');
    %plot(qrs_pos,thI1,'-m');

    hold off;
end