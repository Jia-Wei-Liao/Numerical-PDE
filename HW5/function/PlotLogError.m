function PlotLogError(mList, ErrorList)
loglog(mList, ErrorList, '-ro', 'LineWidth', 1.2);
axis([min(mList), max(mList), ...
    min(ErrorList), max(ErrorList)]);
xlabel('$\log h $', 'interpreter', 'latex');
ylabel('$\log e$', 'interpreter', 'latex');

end