% 获取当前目录下所有.mat文件
files = dir('*.mat');

% 循环处理每个文件
for i = 1:length(files)
    filename = files(i).name;
    
    % 加载.mat文件到结构体
    s = load(filename);
    
    % 检查是否存在'err'字段
    if isfield(s, 'err')
        % 将err中所有0替换为指定值
        s.err(s.err < 9.5e-13) = 1e-12;
        
        % 保存修改后的数据回原文件
        save(filename, '-struct', 's');
    end
end
