% Some big files can't be saved to git. Are they still around?

all_files_are_there = true;

for idx = 1:19
    test_name = strcat("test_setup_assembly.pa", string(idx), ".patxt");
    if exist(test_name, 'file') == 0
        all_files_are_there = false;
        fprintf("*** WARNING ***\n" + ...
                "You are missing some large files that couldn't be saved to git");
        break
    end
end