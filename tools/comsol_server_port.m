function port = comsol_server_port()
value = getenv('COMSOL_SERVER_PORT');
if isempty(value), value = '2036'; end
port = str2double(value);
validateattributes(port, {'double'}, {'scalar','integer','>=',1,'<=',65535});
end
